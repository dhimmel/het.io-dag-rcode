options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'

# Load code
code.dir <- file.path(project.dir, 'rcode')
sources <- c('settings.R', 'hnlp-learning.R', 'machine-learning.R', 'plotting.R')
for (source.filename in sources) {
  source(file.path(code.dir, source.filename))
}

# Load network info
network.id <- '140615-0-training'
network.dir <- file.path(project.dir, 'networks', 'permuted', network.id)
dirs <- InitializeNetworkDir(network.dir)

# Read Features
feature.df <- ReadFeatures(dirs$features)
feature.names <- colnames(feature.df)[-(1:9)]
enhance.df <- read.delim(file.path(dirs$model, 'enhancing-features.txt'))
enhancing.features <- enhance.df[enhance.df$select, 'feature']

# check that no testing positives were included as associations in the network
stopifnot(! any(feature.df$network_status & feature.df$part == 'test'))

# Create training and testing datasets
train.df <- subset(feature.df, part=='train' & (status == 'negative' | network_status))
test.df <- subset(feature.df, part=='test' & status_int != -1)
X.train <- as.matrix(train.df[, feature.names])
X.test <- as.matrix(test.df[, feature.names])
Xe.train <- as.matrix(train.df[, enhancing.features])
Xe.test <- as.matrix(test.df[, enhancing.features])
y.train <- train.df$network_status
y.test <- test.df$status_int

# Fit and test model
fit.ridge <- TrainModel(X=X.train, y=y.train, alpha=0)
test.ridge <- TestModel(cv.model=fit.ridge$cv.model, X=X.test, y=y.test)
SaveFit(fit.ridge, dirs, suffix='-ridge')
SaveTest(test.ridge, dirs, suffix='-ridge')

fit.lasso <- TrainModel(X=X.train, y=y.train, alpha=1)
test.lasso <- TestModel(cv.model=fit.lasso$cv.model, X=X.test, y=y.test)
SaveFit(fit.lasso, dirs, suffix='-lasso')
SaveTest(test.lasso, dirs, suffix='-lasso')

fit.select <- TrainModel(X=Xe.train, y=y.train, alpha=0)
test.select <- TestModel(cv.model=fit.select$cv.model, X=Xe.test, y=y.test)
SaveFit(fit.select, dirs, suffix='-select')
SaveTest(test.select, dirs, suffix='-select')

## Ridge performance
# ROC
roc.df <- rbind(
  cbind(fit.ridge$vtm$roc.df, 'part'='train'),
  cbind(test.ridge$vtm$roc.df, 'part'='test'))
gg.roc <- ggplot(roc.df, aes(fpr, recall, color=part))
gg.roc <- ggROC(gg.roc) + geom_path() +
  scale_color_manual(values=as.character(solarized[c('violet', 'green')]), 
    name='Partition (AUROC)', breaks=c('test', 'train'),
    labels=c(sprintf('Testing (%.3f)', test.ridge$vtm$auroc), 
             sprintf('Training (%.3f)', fit.ridge$vtm$auroc)))

# Testing Precision-Recall Curve
prc.df <- PrunePRC(test.ridge$vtm$threshold.df)
prc.df$panel <- 'Permuted Ridge'
gg.prc <- ggplot(prc.df[nrow(prc.df):1, ], aes(recall, precision, color=threshold))
gg.prc <- ggPRC(gg.prc) +
  annotate('text', x=0.7, y=Inf, 
    label=sprintf('AUPRC == %s', ChrRound(test.ridge$vtm$auprc, 3)), 
    hjust=1, vjust=2, parse=TRUE, size=3.5) +
  facet_grid(panel ~ .)

# Save combined ROC and PRC
path <- file.path(dirs$plots, 'performance-ridge.pdf')
OpenPDF(path, width=width.full, height=2.5)
gridExtra::grid.arrange(gg.roc, gg.prc, nrow=1, widths=c(1, 1.625))
ClosePDF(path)


## Lasso performance
# ROC
roc.df <- rbind(
  cbind(fit.lasso$vtm$roc.df, 'part'='train'),
  cbind(test.lasso$vtm$roc.df, 'part'='test'))
gg.roc <- ggplot(roc.df, aes(fpr, recall, color=part))
gg.roc <- ggROC(gg.roc) + geom_path() +
  scale_color_manual(values=as.character(solarized[c('violet', 'green')]), 
    name='Partition (AUROC)', breaks=c('test', 'train'),
    labels=c(sprintf('Testing (%.3f)', test.lasso$vtm$auroc), 
             sprintf('Training (%.3f)', fit.lasso$vtm$auroc)))

# Testing Precision-Recall Curve
prc.df <- PrunePRC(test.lasso$vtm$threshold.df)
prc.df$panel <- 'Permuted Lasso'
gg.prc <- ggplot(prc.df[nrow(prc.df):1, ], aes(recall, precision, color=threshold))
gg.prc <- ggPRC(gg.prc) +
  annotate('text', x=0.7, y=Inf, 
    label=sprintf('AUPRC == %s', ChrRound(test.lasso$vtm$auprc, 3)), 
    hjust=1, vjust=2, parse=TRUE, size=3.5) +
  facet_grid(panel ~ .)

# Save combined ROC and PRC
path <- file.path(dirs$plots, 'performance-lasso.pdf')
OpenPDF(path, width=width.full, height=2.5)
gridExtra::grid.arrange(gg.roc, gg.prc, nrow=1, widths=c(1, 1.625))
ClosePDF(path)


## Ridge Select performance
# ROC
roc.df <- rbind(
  cbind(fit.select$vtm$roc.df, 'part'='train'),
  cbind(test.select$vtm$roc.df, 'part'='test'))
gg.roc <- ggplot(roc.df, aes(fpr, recall, color=part))
gg.roc <- ggROC(gg.roc) + geom_path() +
  scale_color_manual(values=as.character(solarized[c('violet', 'green')]), 
    name='Partition (AUROC)', breaks=c('test', 'train'),
    labels=c(sprintf('Testing (%.3f)', test.select$vtm$auroc), 
             sprintf('Training (%.3f)', fit.select$vtm$auroc)))

# Testing Precision-Recall Curve
prc.df <- PrunePRC(test.select$vtm$threshold.df)
prc.df$panel <- 'Permuted Ridge Select'
gg.prc <- ggplot(prc.df[nrow(prc.df):1, ], aes(recall, precision, color=threshold))
gg.prc <- ggPRC(gg.prc) +
  annotate('text', x=0.7, y=Inf, 
    label=sprintf('AUPRC == %s', ChrRound(test.select$vtm$auprc, 3)), 
    hjust=1, vjust=2, parse=TRUE, size=3.5) +
  facet_grid(panel ~ .)

# Save combined ROC and PRC
path <- file.path(dirs$plots, 'performance-select.pdf')
OpenPDF(path, width=width.full, height=2.5)
gridExtra::grid.arrange(gg.roc, gg.prc, nrow=1, widths=c(1, 1.625))
ClosePDF(path)


