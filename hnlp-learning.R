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
  test <- list('X'=X, 'y'=y)
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

  pred.df <- data.frame('status'=test$y, 'prediction'=test$y.predicted, test$X, check.names=FALSE)
  path <- file.path(dirs$model, sprintf('testing-predictions%s.txt.gz', suffix))
  gz.file <- gzfile(path, 'w')
  formatted.df <- format(pred.df, digits=digits)
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



# AUCDF Helpers
auc.colnames <- c('breadth', 'type', 'name', 'feature', 
  'metric', 'disease_code', 'disease_name',
  'negatives', 'positives', 'auroc', 'auprc')

OrderAUC <- function(auc.df, blank.df) {
  auc.df <- merge(auc.df, blank.df, all.x=TRUE, all.y=TRUE)
  auc.df <- auc.df[, match(auc.colnames, colnames(auc.df))]
  return(auc.df)
}


GlobalModelAUCDF <- function(fit) {
  # AUCs for a model
  auc.df <- data.frame('breadth'='global', 'type'='model',
    'negatives'=sum(fit$y == 0), 'positives'=sum(fit$y == 1),
    'auroc'=fit$vtm$auroc, 'auprc'=fit$vtm$auprc)
  return(auc.df)
}

SpecificModelAUCDF <- function(fit, disease_codes) {
  # disease specific AUCs for a model
  score.split <- split(fit$y.predicted, disease_codes)
  status.split <- split(fit$y, disease_codes)
  stopifnot(names(score.split) == names(status.split))
  auc.df <- t(mapply(function(score, status) 
    {vtm <- VariableThresholdMetrics(score, status); c(vtm$auroc, vtm$auprc)},
    score=score.split, status=status.split))
  colnames(auc.df) <- c('auroc', 'auprc')
  auc.df <- data.frame(
    'breadth' = 'disease_specific',
    'type'='model',
    'disease_code'=names(score.split),
    'negatives'=sapply(status.split, function(y) sum(y == 0)),
    'positives'=sapply(status.split, function(y) sum(y == 1)),
    auc.df)
  return(auc.df)
}

FeatureAUCDF <- function(X, y) {
  # global AUCs for each feature
  auc.df <- data.frame(
    'feature'=colnames(X),
    'negatives'=sum(y == 0),
    'positives'=sum(y == 1))
  auc.df <- cbind(auc.df, 
    t(apply(X, 2, function(x) {
      vtm <- VariableThresholdMetrics(x, y)
      return(c(vtm$auroc, vtm$auprc))})) )
  colnames(auc.df) <- c('feature', 'negatives', 'positives', 'auroc', 'auprc')
  return(auc.df)
}

GlobalFeatureAUCDF <- function(X, y) {
  # global AUCs for each feature
  auc.df <- data.frame(
    FeatureAUCDF(X, y),
    'type'='feature',
    'breadth'='global')
  return(auc.df)
}

HelperSpecificFeatureAUCDF <- function(Xy.df, feature.names) {
  auc.df <- FeatureAUCDF(X=Xy.df[, feature.names], y=Xy.df$status)
  return(auc.df)
}

SpecificFeatureAUCDF <- function(X, y, disease_codes) {
  # disease-specific AUCs for each feature
  feature.names <- colnames(X)
  Xy.df <- data.frame('disease_code'=disease_codes, 'status'=y, X, check.names=FALSE)
  auc.df <- data.frame(
    plyr::ddply(Xy.df, 'disease_code', HelperSpecificFeatureAUCDF, feature.names=feature.names),
    'type'='feature',
    'breadth'='disease_specific')
  return(auc.df)
}


ComputeAUCDF <- function(X, y, disease_codes, fit.list) {

  # Create blank data.frame
  blank.df <- data.frame(t(rep(NA, length(auc.colnames))))
  names(blank.df) <- auc.colnames
  blank.df <- blank.df[-1, ]

  # Calculate global model AUCs
  global.model.df <- plyr::ldply(fit.list, GlobalModelAUCDF, .id='name')
  global.model.df <- OrderAUC(global.model.df, blank.df)
  
  # Calculate disease-specific model AUCs
  specific.model.df <- plyr::ldply(fit.list, SpecificModelAUCDF, disease_codes=disease_codes, .id='name')
  specific.model.df <- OrderAUC(specific.model.df, blank.df)
  
  # Calculate global feature AUCs
  global.feature.df <- GlobalFeatureAUCDF(X, y)
  global.feature.df <- OrderAUC(global.feature.df, blank.df)

  # Calculate disease-specific feature AUCs
  specific.feature.df <- SpecificFeatureAUCDF(X, y, disease_codes)
  specific.feature.df <- OrderAUC(specific.feature.df, blank.df)

  # Combine and convert factor columns to character
  auc.df <- rbind(global.model.df, specific.model.df, global.feature.df, specific.feature.df)
  factor.col <- sapply(auc.df, is.factor)
  auc.df[factor.col] <- lapply(auc.df[factor.col], as.character)

  # Add columns with additional information
  auc.df$disease_pathophys <- pathophys.df[
    match(auc.df$disease_code, pathophys.df$disease_code), 'pathophysiology']
  auc.df$disease_name <- pathophys.df[
    match(auc.df$disease_code, pathophys.df$disease_code), 'disease_name']
  auc.df$metric <- desc.df[match(auc.df$feature, desc.df$feature), 'metric']
  is.feature <- auc.df$type == 'feature'
  auc.df[is.feature, 'name'] <- feature.converter[auc.df[is.feature, 'feature']]

  return(auc.df)
}

