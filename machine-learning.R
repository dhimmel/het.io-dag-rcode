SplitTesting <- function(tab, status, percent.testing=1/3, seed=0) {
  # Split dataframe or matrix tab with binary status into train and test sets.
  # Dataset is stratified by status for equal percent.testing in both sets.
  set.seed(seed)
  status <- as.logical(status)
  tab.pos <- tab[status, ]
  tab.neg <- tab[! status, ]
  number.pos <- nrow(tab.pos)
  number.neg <- nrow(tab.neg)
  number.testing.pos <- round(number.pos * percent.testing)
  number.testing.neg <- round(number.neg * percent.testing)
  number.training.pos <- number.pos - number.testing.pos
  number.training.neg <- number.neg - number.testing.neg

  pos.is.testing <- rep(c(TRUE, FALSE), c(number.testing.pos, number.training.pos))
  neg.is.testing <- rep(c(TRUE, FALSE), c(number.testing.neg, number.training.neg))
  pos.is.testing <- sample(pos.is.testing)
  neg.is.testing <- sample(neg.is.testing)
  testing <- rbind(tab.pos[pos.is.testing, ], tab.neg[neg.is.testing, ])
  training <- rbind(tab.pos[! pos.is.testing, ], tab.neg[! neg.is.testing, ])
  tab.list <- list('training'=training, 'testing'=testing)
  return(tab.list)
}

VariableThresholdMetrics <- function(score, status) {
  # TPR is equivalent to recall
  rocr.pred <- ROCR::prediction(score, status)
  auroc <- ROCR::performance(rocr.pred, 'auc')@y.values[[1]]
  threshold.df <- data.frame(
    'threshold'=rocr.pred@cutoffs[[1]],
    'fpr'=ROCR::performance(rocr.pred, measure='fpr')@y.values[[1]],
    'recall'=ROCR::performance(rocr.pred, measure='rec')@y.values[[1]],
    'precision'=ROCR::performance(rocr.pred, measure='prec')@y.values[[1]],
    'lift'=ROCR::performance(rocr.pred, measure='lift')@y.values[[1]]
  )

  roc.df <- threshold.df[, c('fpr', 'recall')]
  for (measure in c('fpr', 'recall')) {
    not.dup <- ! duplicated(roc.df$recall)
    not.dup <- not.dup | c(not.dup[-1], TRUE)
    roc.df <- roc.df[not.dup, ]
  }

  metrics <- list('auroc'=auroc, 'threshold.df'=threshold.df, 'roc.df'=roc.df)
  return(metrics)
}

PrunePRC <- function(prc.df, min.dist=0.0005) {
  dist.df <- prc.df[, c('precision', 'recall')]
  pointer <- 1
  as.index <- sapply(2:nrow(dist.df), function(i) {
    distance <- dist(dist.df[c(pointer, i), 1:2])[1]
    if (distance > min.dist) {
      pointer <<- i
      return(i)
    } else {return(pointer)}
  })
  prc.df <- prc.df[c(1, unique(as.index)), ]
  return(prc.df)
}

CrossValidationFolds <- function(Y, nfolds=10, seed=0) {
  stopifnot(min(table(Y)) >= nfolds)
  positions <- seq(1, length(Y))
  set.seed(seed)
  neg.positions <- sample(positions[Y == 0])
  pos.positions <- sample(positions[Y == 1])
  neg.fold.assigs <- rep_len(1:nfolds, length.out=length(neg.positions))
  pos.fold.assigs <- rep_len(1:nfolds, length.out=length(pos.positions))
  neg.folds <- split(neg.positions, neg.fold.assigs)
  pos.folds <- split(pos.positions, pos.fold.assigs)
  folds <- mapply(c, neg.folds, pos.folds, SIMPLIFY=FALSE)
  return(folds)
}

