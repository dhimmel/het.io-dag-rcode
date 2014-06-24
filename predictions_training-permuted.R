library(glmnet)
library(ggplot2)
library(doMC)
library(grid)
library(gtools)
library(plyr)
library(reshape2)
library(gridExtra)

options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
network.id <- '140615-training-permuted'

code.dir <- file.path(project.dir, 'rcode')
network.dir <- file.path(project.dir, 'networks', network.id)
feature.dir <- file.path(network.dir, 'features')
graphics.dir <- file.path(network.dir, 'graphics')
modeling.dir <- file.path(network.dir, 'modeling')
dir.create(graphics.dir, showWarnings=FALSE)
dir.create(modeling.dir, showWarnings=FALSE)

source(file.path(code.dir, 'machine-learning.R'))
source(file.path(code.dir, 'functions.R'))

glmnet.alpha <- 0 # 0 for ridge regression

# Read Features
feature.df <- ReadFeatures(feature.dir)
feature.names <- colnames(feature.df)[-(1:9)]

stopifnot(! any(feature.df$network_status & feature.df$part == 'test'))

train.df <- subset(feature.df, part=='train' & (status == 'negative' | network_status))
test.df <- subset(feature.df, part=='test' & status_int != -1)

X.train <- as.matrix(train.df[, feature.names])
X.test <- as.matrix(test.df[, feature.names])
y.train <- train.df$network_status
y.test <- test.df$status_int

doMC::registerDoMC(cores=7)
set.seed(0)
cv.ridge.train <- glmnet::cv.glmnet(X.train, y.train, family='binomial', 
  alpha=glmnet.alpha, standardize=TRUE, parallel=TRUE)
# Save cv.glmnet object as an R data structure
saveRDS(cv.ridge.train, file.path(modeling.dir, 'training-model.rds'))
#cv.ridge.train <- readRDS(file.path(modeling.dir, 'training-model.rds'))
lambda.train <- cv.ridge.train$lambda.1se
y.predicted.train <- as.numeric(predict(cv.ridge.train, s=lambda.train, newx=X.train, type='response'))
y.predicted.test <- as.numeric(predict(cv.ridge.train, s=lambda.train, newx=X.test, type='response'))

vtm.train <- VariableThresholdMetrics(y.predicted.train, y.train)
vtm.test <- VariableThresholdMetrics(y.predicted.test, y.test)
vtm.test.path <- file.path(modeling.dir, 'variable-threshold-metrics-testing.txt')
write.table(round(vtm.test$threshold.df, 5), vtm.test.path, sep='\t', row.names=FALSE, quote=FALSE)

# Ceofficient data.frame
feature.name.mat <- do.call(rbind, strsplit(feature.names, '[|]'))
coef.df <- GLMNetCoef(cv.ridge.train, X.train, y.train, 
  metric=c(NA, feature.name.mat[, 1]), metapath=c(NA, feature.name.mat[, 2]))
coef.path <- file.path(modeling.dir, 'coefficients.txt')
write.table(coef.df, coef.path, sep='\t', row.names=FALSE, quote=FALSE)

# Testing Precision-Recall Curve
prc.df <- PrunePRC(vtm.test$threshold.df)
gg.prc <- ggplot(prc.df[nrow(prc.df):1, ], aes(recall, precision, color=threshold))
gg.prc <- ggPRC(gg.prc)

# Testing ROC Curve
roc.df <- rbind(cbind(vtm.train$roc.df, 'part'='train'), cbind(vtm.test$roc.df, 'part'='test'))
write.table(roc.df, file.path(modeling.dir, 'ROC-testing-and-training.txt'), sep='\t', row.names=FALSE, quote=FALSE)

gg.roc <- ggplot(roc.df, aes(fpr, recall, color=part))
gg.roc <- ggROC(gg.roc) + geom_path() +
  scale_color_manual(values=as.character(solarized[c('violet', 'green')]), 
    name='Partition (AUROC)', breaks=c('test', 'train'),
    labels=c(sprintf('Testing (%.3f)', vtm.test$auroc), 
             sprintf('Training (%.3f)', vtm.train$auroc)))

# Save combined ROC and PRC
path <- file.path(graphics.dir, 'performance-permutation-testing.pdf')
OpenPDF(path, width=width.full, height=2.5)
gridExtra::grid.arrange(gg.roc, gg.prc, nrow=1, widths=c(1, 1.625))
ClosePDF(path)



