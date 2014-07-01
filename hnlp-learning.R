library(doMC)
library(glmnet)


ReadFeatures <- function(feature.dir) {
  # Read and combine the feature files in feature.dir with names beginning with 'DOID_'.
  feature.filenames <- list.files(feature.dir, pattern='DOID_')
  feature.df.list <- list()
  for (feature.filename in feature.filenames) {
    doid_code <- gsub('.txt.gz', '', feature.filename)
    cat(sprintf('Reading features for %s\n', doid_code))
    feature.path <- file.path(feature.dir, feature.filename)
    feature.df <- read.delim(feature.path, check.names=FALSE)
    feature.df.list[[doid_code]] <- feature.df
  }
  return(do.call(rbind, feature.df.list))
}


GLMNetCoef <- function(cv.fit, X, y, prepend='', ...) {
  lambda <- cv.fit$lambda.1se
  coef.vec <- coef(cv.fit, s=lambda)[, 1]
  zcoef.vec <- c(0, coef.vec[-1] * apply(X, 2, sd) / sd(y))
  coef.df <- data.frame('feature'=c('intercept', colnames(X)), ...,
    'coef'=coef.vec, 'zcoef'=zcoef.vec, row.names=NULL)
  colnum <- ncol(coef.df)
  colnames(coef.df)[c(colnum - 1, colnum)] <- paste0(prepend, c('coef', 'zcoef'))
  return(coef.df)
}


TrainModel <- function(X, y, alpha, cores=7, seed=0) {
  # Fit a regularized logistic regression model using the glmnet package.
  # alpha is the regularization parameter (0 for ridge, 1 for lasso).
  fit <- list('X'=X, 'y'=y, 'alpha'=alpha, 'seed'=seed)
  
  # train model
  doMC::registerDoMC(cores=cores)
  set.seed(seed)
  fit$cv.model <- cv.fit <- glmnet::cv.glmnet(X, y, family='binomial', 
    alpha=alpha, standardize=TRUE, parallel=TRUE)
  fit$lambda.1se <- lambda.1se <- cv.fit$lambda.1se

  # model information and performance
  fit$coef.df <- GLMNetCoef(cv.fit, X, y)
  fit$y.predicted <- y.predicted <- as.numeric(predict(cv.fit, s=lambda.1se, newx=X, type='response'))
  fit$vtm <- VariableThresholdMetrics(y.predicted, y)
  
  return(fit)
}

TestModel <- function(cv.model, X, y) {
  test <- list()
  test$y.predicted <- y.predicted <- as.numeric(
    predict(cv.model, s=cv.model$lambda.1se, newx=X, type='response'))
  test$vtm <- vtm <- VariableThresholdMetrics(y.predicted, y)
  return(test)
}

MakePredictions <- function(cv.model, X) {
  y.predicted <- as.numeric(
    predict(cv.model, s=cv.model$lambda.1se, newx=X, type='response'))
  return(y.predicted)
}

SaveFit <- function(fit, dirs, suffix='', digits=5) {
  # Save cv.model
  path <- file.path(dirs$model, sprintf('model%s.rds', suffix))
  saveRDS(fit$cv.model, path)

  # Save training vtm
  path <- file.path(dirs$model, sprintf('training-vtm%s.txt.gz', suffix))
  gz.file <- gzfile(path, 'w')
  formatted.df <- format(fit$vtm$threshold.df, digits=digits)
  write.table(formatted.df, gz.file, sep='\t', row.names=FALSE, quote=FALSE)
  close(gz.file)

  # Save coef.df
  path <- file.path(dirs$model, sprintf('coefficients%s.txt', suffix))
  formatted.df <- format(fit$coef.df, digits=digits)
  write.table(formatted.df, path, sep='\t', row.names=FALSE, quote=FALSE)
}

SaveTest <- function(test, dirs, suffix='', digits=5) {
  # Save testing vtm
  path <- file.path(dirs$model, sprintf('testing-vtm%s.txt.gz', suffix))
  gz.file <- gzfile(path, 'w')
  formatted.df <- format(test$vtm$threshold.df, digits=digits)
  write.table(formatted.df, gz.file, sep='\t', row.names=FALSE, quote=FALSE)
  close(gz.file)
}

InitializeNetworkDir <- function(network.dir) {
  dirs <- list()
  # create a list of directory paths
  for (dirname in c('features', 'plots', 'model')) {
    dirs[[dirname]] <- file.path(network.dir, dirname)
  }
  # make directories if non-existent
  for (dirname in c('plots', 'model')) {
    dir.create(dirs[[dirname]], showWarnings=FALSE)
  }
  return(dirs)
}
