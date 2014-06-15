IHS <- function(x, theta){
  # Inverse Hyperbolic Sine Transformation
  # http://goo.gl/I4pDoq
  x_transformed <- (1 / theta) * asinh(theta * x)
  return(x_transformed)
}

LogLikIHS <- function(theta, x){
  # IHS transform x and return the log likelihood of the resulting distribution
  # http://goo.gl/I4pDoq
  n <- length(x)
  x_transformed <- IHS(x, theta)
  x_transformed_mean <- mean(x_transformed)
  x_transformed_totvar <- sum((x_transformed - x_transformed_mean)^2)
  log.lik <- -n * log(x_transformed_totvar) - sum(log(1 + theta^2 * x^2))
  return(log.lik)
}


FitThetaIHS <- function(x, lower=0.001, upper=50) {
  # Return optimized thetas for an IHS transformation for normality and likelihood methods.
  #shapiro.opt <- optimize(ShapiroIHS, lower=lower, upper=upper, x=x, maximum=TRUE)
  loglik.opt <- optimize(LogLikIHS, lower=lower, upper=upper, x=x, maximum=TRUE)
  #ks.opt <- optimize(ksIHS, lower=lower, upper=upper, x=x, maximum=TRUE)
  #thetas <- c('shapiro'=shapiro.opt$maximum, 'ks'=ks.opt, 'loglik'=loglik.opt$maximum)
  thetas <- c('loglik'=loglik.opt$maximum)
  return(thetas)
}


library(glmnet)
library(ggplot2)
library(car)
library(reshape2)
options(stringsAsFactors=FALSE)
source('machine-learning.R')


project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
feature.dir <- file.path(project.dir, 'networks', 'elderbatch', '140514-all-assoc', 'features')
feature.filenames <- list.files(feature.dir, pattern='DOID_')
feature.df.list <- list()
for (feature.filename in feature.filenames) {
  doid_code <- gsub('.txt.gz', '', feature.filename)
  cat(sprintf('Reading features for %s\n', doid_code))
  feature.path <- file.path(feature.dir, feature.filename)
  feature.df <- read.delim(feature.path, check.names=FALSE)
  feature.df.list[[doid_code]] <- feature.df
}
feature.df <- do.call(rbind, feature.df.list)
feature.names <- colnames(feature.df)[-(1:8)]


train.feature.df <- subset(feature.df, status_int >= 0)
X <- as.matrix(train.feature.df[, feature.names])
y <- train.feature.df$status_int

X.ihs <- apply(X, 2, function(x) IHS(x, FitThetaIHS(x)['loglik']))
X.comb <- cbind(X[, 1:2], X.ihs[, -(1:2)])

X.list <- list('untransformed'=X, 'ihs'=X.ihs, 'combined'=X.comb)
auroc.list <- list()
prediction.list <- list()

doMC::registerDoMC(cores=7)
glmnet.alpha <- 0 # ridge regression
for (trans.name in names(X.list)) {
  cat(sprintf('Beginning Work on %s\n', trans.name))
  XX <- X.list[[trans.name]]
  set.seed(0)
  cv.ridge <- glmnet::cv.glmnet(XX, y, family='binomial', 
    alpha=glmnet.alpha, standardize=TRUE, parallel=TRUE)
  lambda.global <- cv.ridge$lambda.1se
  y.predicted <- as.numeric(predict(cv.ridge, s=lambda.global, newx=XX, type='response'))
  auroc.list[[trans.name]] <- VariableThresholdMetrics(y.predicted, y)$auroc
  prediction.list[[trans.name]] <- y.predicted
}





prediction.df <- as.data.frame(do.call(cbind, prediction.list))
#prediction.melt <- reshape2::melt(prediction.df, variable.name='transformation', value.name='prediction')
#ggplot(prediction.melt, aes(transformation, prediction)) +
  



x.melt <- data.frame(reshape2::melt(as.data.frame(X), variable.name='feature'), 'transformation'='none')
x.ihs.melt <- data.frame(reshape2::melt(as.data.frame(X.ihs), variable.name='feature'), 'transformation'='ihs')
melt.df <- rbind(x.melt, x.ihs.melt)


split.df <- split(melt.df, melt.df$feature)
pdf('~/Desktop/gghist.pdf', width=5, height=7)
for (partial.df in split.df) {
  gghist <- ggplot(partial.df, aes(value)) + geom_histogram() + 
    facet_wrap( ~ transformation, scale='free', ncol=1) + 
    scale_y_sqrt() + theme_bw() + ggtitle(partial.df[1, 'feature'])
  print(gghist)
}
dev.off()
#ggsave('~/Desktop/gghist.pdf', gghist, width=8, height=100, limitsize=FALSE)

