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
directory <- file.path(project.dir, 'networks', '140518-metricsweep')

source(file.path(code.dir, 'machine-learning.R'))
# Using python processing
#aucs.path <- file.path(directory, 'feature-aucs.txt')
#auc.df <- read.delim(aucs.path)
#auc.df$rank <- rank(auc.df$auc)
#ggplot(auc.df, aes(metapath, metric, fill=rank)) + geom_tile()
#ggplot(auc.df, aes(metric, auc)) + geom_point()

# Using R processing
feature.path <- file.path(directory, 'features.txt.gz')
feature.df <- read.delim(feature.path, check.names=FALSE, stringsAsFactors=FALSE)

# Read partition file with master list of testing and training
# Not needed because non-training pairs are not computed, left in for safety.
partitions.path <- file.path(project.dir, 'partitions.txt.gz')
part.df <- read.delim(partitions.path, colClasses='character')
# Filter feature.df to only include training pairs
training.pairs <- apply(part.df[part.df$part=='train', c('disease_code', 'gene_symbol')], 1, paste, collapse='|')
feature.df.pairs <- apply(feature.df[, c('target', 'source')], 1, paste, collapse='|')
feature.df <- feature.df[feature.df.pairs %in% training.pairs, ]

# Additional feature.df processing
feature.df[is.na(feature.df)] <- 0
feature.names <- colnames(feature.df)[-seq(4)]
stat.df <- as.data.frame(do.call(rbind, strsplit(feature.names, ':')))
colnames(stat.df) <- c('metric', 'metapath')
stat.df$feature <- feature.names
stat.df <- subset(stat.df, metapath != 'G-a-D')
#metrics <- setdiff(stat.df$metric, c('PCs', 'PCt'))
metrics <- setdiff(stat.df$metric, c('PC', 'PCs', 'PCt', 'DWPC_0.9', 'DWPC_1.0'))

doMC::registerDoMC(cores=8)

glmnet.alpha <- 0.0
kfolds <- 20
repititions <- 10
metric.df <- data.frame('metric'=metrics,
  'fold'=rep(1:kfolds, each=length(metrics)),
  'repitition'=rep(1:repititions, each=length(metrics)*kfolds),
  'lambda'=NA, 'auc'=NA)
for (i in 1:repititions) {
  folds <- CrossValidationFolds(feature.df$status, kfolds, seed=i)

  for (k in 1:kfolds) {
    cat(sprintf('Repition %s, Fold %s\n', i, k))
    fold <- folds[[k]]

    y.train <- feature.df$status[-fold]
    y.test <- feature.df$status[fold]

    for (metric in metrics) {
      cat(sprintf('Metric: %s\n', metric))
      features <- stat.df[stat.df$metric == metric, 'feature']
      features <- c(features, 'PCs:G-a-D', 'PCt:G-a-D')
      X.train <- as.matrix(feature.df[-fold, features])
      X.test <- as.matrix(feature.df[fold, features])
      cv.enet <- glmnet::cv.glmnet(X.train, y.train, family='binomial', alpha=glmnet.alpha, standardize=TRUE, parallel=TRUE)
      lambda <- cv.enet$lambda.1se
      coef(cv.enet, lambda)
      y.predicted <- predict(cv.enet, s=lambda, newx=X.test, type='response')
      vtm <- VariableThresholdMetrics(y.predicted, y.test)
      metric.df.index <- metric.df$metric == metric & metric.df$fold == k & metric.df$repitition == i
      metric.df[metric.df.index, 'lambda'] <- lambda
      metric.df[metric.df.index, 'auc'] <- vtm$auroc
      # Precision-recall Curve
      #ggplot(vtm$prg, aes(recall, precision, color=threshold)) + 
      #  geom_line(size=2) + theme_bw() + scale_colour_gradientn(colours = rainbow(100))
    }
  }
}




enet.auc.path <- file.path(directory, sprintf('ridge-auc-%sx%s-CV.txt', repititions, kfolds))
write.table(metric.df, enet.auc.path, sep='\t', row.names=FALSE, quote=FALSE)
#confint on ggplot
#mean.conf.int <- function(x) {t.test(x)$conf.int[1:2]}
#stat_summary(fun.y='mean.conf.int', geom='line')

metric.df <- read.delim(enet.auc.path, stringsAsFactors=FALSE)
metric.summary.df <- plyr::ddply(metric.df, c('metric', 'repitition'),
  plyr::summarize, 'auc'=mean(auc))

metric.summary.df[, 'dwpc_exponent'] <- suppressWarnings(as.numeric(gsub('DWPC_', '', metric.summary.df$metric)))
npc.aucs <- subset(metric.summary.df, metric == 'NPC')[, 'auc']
npc.mean <- mean(npc.aucs)
npc.confint <- t.test(npc.aucs)$conf.int

threshold.dwpc <- 0.4

dwpc.plot.df <- na.omit(metric.summary.df)
ggplot(dwpc.plot.df, aes(dwpc_exponent, auc)) + theme_bw() +
  geom_vline(xintercept=threshold.dwpc, color='darkgreen', linetype='dashed') + 
  geom_rect(xmin=-Inf, xmax=Inf, ymin=npc.confint[1], ymax=npc.confint[2], fill='#7EC0EE') +
  geom_hline(yintercept=npc.mean, color='#42647F') +
  #geom_point() +
  geom_smooth(linetype=0, method='loess', span=0.75, alpha=0.8) +
  stat_summary(fun.y='mean', geom='point') + 
  #theme(axis.text.x = element_text(angle=60, hjust=1)) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) +
  #coord_cartesian(ylim=c(0.784, 0.821)) +
  xlab('DWPC Exponent') + ylab(expression(10 %*% 20 ~ fold ~ CV ~ AUROC))


enet.auc.plot.path <- file.path(directory, sprintf('ridge-auc-%sx%s-CV.pdf', repititions, kfolds))
ggsave(enet.auc.plot.path, width=3, height=3)



