options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'

# Load code
code.dir <- file.path(project.dir, 'rcode')
sources <- c('settings.R', 'hnlp-learning.R', 'machine-learning.R', 'plotting.R')
for (source.filename in sources) {
  source(file.path(code.dir, source.filename))
}

# Load network info
network.id <- '140615-training'
network.dir <- file.path(project.dir, 'networks', network.id)
dirs <- InitializeNetworkDir(network.dir)

# Read Features
feature.df <- ReadFeatures(dirs$features)
feature.names <- colnames(feature.df)[-(1:8)]

# Create training and testing datasets
train.df <- subset(feature.df, part=='train' & status_int != -1)
test.df <- subset(feature.df, part=='test' & status_int != -1)
X.train <- as.matrix(train.df[, feature.names])
X.test <- as.matrix(test.df[, feature.names])
y.train <- train.df$status_int
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


# Training and Testing ROC Curves

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
prc.df$panel <- 'Ridge'
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
  scale_color_manual(values=Solar('violet', 'green')), 
    name='Partition (AUROC)', breaks=c('test', 'train'),
    labels=c(sprintf('Testing (%.3f)', test.lasso$vtm$auroc), 
             sprintf('Training (%.3f)', fit.lasso$vtm$auroc)))

# Testing Precision-Recall Curve
prc.df <- PrunePRC(test.lasso$vtm$threshold.df)
prc.df$panel <- 'Lasso'
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




