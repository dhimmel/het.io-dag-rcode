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



# AUROCDF Helpers
auroc.colnames <- c('breadth', 'type', 'name', 'feature', 
  'metric', 'disease_code', 'disease_name',
  'positives', 'auroc')

OrderAurocDf <- function(auroc.df, blank.df) {
  auroc.df <- merge(auroc.df, blank.df, all.x=TRUE, all.y=TRUE)
  auroc.df <- auroc.df[, match(auroc.colnames, colnames(auroc.df))]
  return(auroc.df)
}

ComputeFeatureAUC <- function(feat.df) {
  apply(feat.df[, feature.names], 2, function(x) VariableThresholdMetrics(x, feat.df$status_int)$auroc)
}


ComputeAUROCDF <- function(feat.perf.df) {

  # Create blank data.frame
  blank.df <- data.frame(t(rep(NA, length(auroc.colnames))))
  names(blank.df) <- auroc.colnames
  blank.df <- blank.df[-1, ]

  # Calculate disease specific AUCs from the global ridge and lasso models
  auroc.df1 <- rbind(
  plyr::ddply(feat.perf.df, c('disease_code'), plyr::summarize,
    'breadth' = 'disease_specific',
    'type'='model', 'name'='Ridge',
    'auroc' = VariableThresholdMetrics(ridge, status_int)$auroc,
    'positives'=sum(status_int)),
  plyr::ddply(feat.perf.df, c('disease_code'), plyr::summarize,
    'breadth' = 'disease_specific',
    'type'='model', 'name'='Lasso',
    'auroc' = VariableThresholdMetrics(lasso, status_int)$auroc,
    'positives'=sum(status_int))
  )
  auroc.df1 <- OrderAurocDf(auroc.df1, blank.df)

  # Calculate global AUCs for each feature
  auroc.vec2 <- ComputeFeatureAUC(feat.perf.df)
  auroc.df2 <- data.frame(
    'feature'=names(auroc.vec2), 
    'auroc'=as.numeric(auroc.vec2),
    'type'='feature',
    'breadth'='global',
    'positives'=sum(feat.perf.df$status_int)
  )
  auroc.df2 <- OrderAurocDf(auroc.df2, blank.df)

  # Calculate disease-specific AUROCs for each feature
  fXd.auroc.tab <- plyr::ddply(feat.perf.df, c('disease_code'), ComputeFeatureAUC)
  auroc.df3 <- reshape2::melt(fXd.auroc.tab, id.vars=c('disease_code'), variable.name='feature', value.name='auroc')
  auroc.df3 <- cbind(auroc.df3, 'type'='feature', 'breadth'='disease_specific')


  # Global AUC for Ridge Model
  auroc.df4 <- rbind(
  data.frame('breadth'='global', 'type'='model', 'name'='Ridge',
    'positives'=sum(feat.perf.df$status_int), 'auroc'=vtm.ridge$auroc),
  data.frame('breadth'='global', 'type'='model', 'name'='Lasso',
    'positives'=sum(feat.perf.df$status_int), 'auroc'=vtm.lasso$auroc)
  )
  auroc.df4 <- OrderAurocDf(auroc.df4, blank.df)

  # Add positives to auroc.df3
  auroc.df3$positives <- auroc.df1[match(auroc.df3$disease_code, auroc.df1$disease_code), 'positives']
  auroc.df3 <- OrderAurocDf(auroc.df3, blank.df)

  # Combine
  auroc.df <- rbind(auroc.df1, auroc.df2, auroc.df3, auroc.df4)
  auroc.df$disease_pathophys <- pathophys.df[match(auroc.df$disease_code, pathophys.df$disease_code), 'pathophysiology']
  auroc.df$disease_name <- pathophys.df[match(auroc.df$disease_code, pathophys.df$disease_code), 'disease_name']

  auroc.df$metric <- desc.df[match(auroc.df$feature, desc.df$feature), 'metric']
  auroc.df[is.na(auroc.df$name), 'name'] <- feature.converter[auroc.df[is.na(auroc.df$name), 'feature']]

  return(auroc.df)
}

