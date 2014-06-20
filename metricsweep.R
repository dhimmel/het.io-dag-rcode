library(ggplot2)
library(gridExtra)
library(glmnet)
library(doMC)
library(ROCR)
library(grid)
library(plyr)

options(stringsAsFactors=FALSE)


project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
code.dir <- file.path(project.dir, 'rcode')
directory <- file.path(project.dir, 'networks', '140615-metricsweep')

source(file.path(code.dir, 'machine-learning.R'))
source(file.path(code.dir, 'functions.R'))

# Using R processing
feature.path <- file.path(directory, 'features.txt.gz')
feature.df <- read.delim(feature.path, check.names=FALSE, stringsAsFactors=FALSE)

# Additional feature.df processing
feature.df[is.na(feature.df)] <- 0
feature.names <- colnames(feature.df)[-seq(4)]
stat.df <- as.data.frame(do.call(rbind, strsplit(feature.names, ':')))
colnames(stat.df) <- c('metric', 'metapath')
stat.df$feature <- feature.names
stat.df <- subset(stat.df, metapath != 'G-a-D')
#metrics <- setdiff(stat.df$metric, c('PCs', 'PCt'))
metrics <- setdiff(stat.df$metric, c('PC', 'PCs', 'PCt', 'DWPC_0.9', 'DWPC_1.0'))

doMC::registerDoMC(cores=7)


glmnet.alphas <- seq(0, 1, 0.2)
kfolds <- 20
repititions <- 10

metric.df.list <- list()
for (glmnet.alpha in glmnet.alphas) {

  for (i in 1:repititions) {
    folds <- CrossValidationFolds(feature.df$status, kfolds, seed=i)

    for (k in 1:kfolds) {
      cat(sprintf('Repition %s, Fold %s\n', i, k))
      fold <- folds[[k]]

      y.train <- feature.df$status[-fold]
      y.test <- feature.df$status[fold]

      set.seed(i)
      for (metric in metrics) {
        cat(sprintf('Metric: %s\n', metric))
        features <- stat.df[stat.df$metric == metric, 'feature']
        features <- c(features, 'PCs:G-a-D', 'PCt:G-a-D')
        X.train <- as.matrix(feature.df[-fold, features])
        X.test <- as.matrix(feature.df[fold, features])
        cv.enet <- glmnet::cv.glmnet(X.train, y.train, family='binomial', 
          alpha=glmnet.alpha, standardize=TRUE, parallel=TRUE)
        lambda <- cv.enet$lambda.1se
        y.predicted <- predict(cv.enet, s=lambda, newx=X.test, type='response')
        vtm <- VariableThresholdMetrics(y.predicted, y.test)
        temp.df <- data.frame('metric'=metric, 'alpha'=glmnet.alpha, 'fold'=k,
          'repitition'=i, 'lambda'=lambda, 'auroc'=vtm$auroc)
        metric.df.list[[length(metric.df.list) + 1]] <- temp.df
      }
    }
  }
}

metric.df <- do.call(rbind, metric.df.list)


metric.path <- file.path(directory, sprintf('model-auc-%sx%s-CV.txt', repititions, kfolds))
write.table(metric.df, metric.path, sep='\t', row.names=FALSE, quote=FALSE)

metric.df <- read.delim(metric.path, stringsAsFactors=FALSE)

metric.summary.df <- plyr::ddply(metric.df, c('alpha', 'metric', 'repitition'),
  plyr::summarize, 'auroc'=mean(auroc))
metric.summary.df[, 'dwpc_exponent'] <- suppressWarnings(as.numeric(gsub('DWPC_', '', metric.summary.df$metric)))

npc.ggdf <- subset(metric.summary.df, metric == 'NPC')
npc.ggdf <- plyr::ddply(npc.ggdf, 'alpha', plyr::summarize,
  'lower'=t.test(auroc, conf.level=0.9999)$conf.int[1],
  'upper'=t.test(auroc, conf.level=0.9999)$conf.int[2],
  'auroc'=mean(auroc))


threshold.dwpc <- 0.4

dwpc.ggdf <- subset(metric.summary.df, ! is.na(dwpc_exponent))
plot.alphas <- unique(dwpc.ggdf$alpha)
vline.df <- data.frame('threshold'=c(threshold.dwpc, rep(NA, length(plot.alphas) - 1)), 'alpha'=plot.alphas)

gg.metric <- ggplot(dwpc.ggdf, aes(dwpc_exponent, auroc)) 
gg.metric <- SetGGTheme(gg.metric) +
  facet_wrap(~ alpha, nrow=3) +
  geom_vline(data=vline.df, aes(xintercept=threshold),
    color=Solar('yellow'), linetype='dashed') + 
  geom_rect(data=npc.ggdf, aes(x=NULL, ymin=lower, ymax=upper), xmin=-Inf, xmax=Inf,
    fill=Solar('violet'), alpha=0.4) +
  geom_hline(data=npc.ggdf, aes(x=NULL, yintercept=auroc), color=Solar('violet'), size=0.25) +
  geom_smooth(linetype=0, method='loess', span=0.75, alpha=0.6,
    fill=Solar('magenta'), level=0.9999) +
  stat_summary(fun.y='mean', geom='point', color=Solar('magenta')) + 
  scale_x_continuous(breaks=seq(0.1, 0.7, 0.2)) +
  scale_y_continuous(breaks=seq(0.76, 0.80, 0.02)) + 
  #coord_cartesian(ylim=c(0.784, 0.821)) +
  xlab('DWPC Exponent') + ylab(expression(10 %*% 20 ~ fold ~ CV ~ AUROC))
#gg.metric

ExpressionFromAlpha <- function(a) {
  if (a == 0) {
    method.name <- '(ridge)'
  } else if (a == 1) {
    method.name <- '(lasso)'
  } else {
    method.name <- '(elastic net)'
  }
  expr.call <- bquote(alpha == .(a) ~ .(method.name))
  return(expr.call)
}

label.list <- lapply(unique(dwpc.ggdf$alpha), ExpressionFromAlpha)
gg.metric <- FacetWrapLabeller(gg.metric, label.list)


enet.auc.plot.path <- file.path(directory, sprintf('metrics-%sx%s-CV.pdf', repititions, kfolds))
ggsave(enet.auc.plot.path, gg.metric, width=width.half, height=3.6)

