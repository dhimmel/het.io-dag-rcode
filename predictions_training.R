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
network.id <- '140615-training'

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
feature.names <- colnames(feature.df)[-(1:8)]

train.df <- subset(feature.df, part=='train' & status_int != -1)
test.df <- subset(feature.df, part=='test' & status_int != -1)

X.train <- as.matrix(train.df[, feature.names])
X.test <- as.matrix(test.df[, feature.names])
y.train <- train.df$status_int
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

coefs.training <- coef(cv.ridge.train, s=lambda.train)[, 1]
zcoefs.training <- c(0, coefs.training[-1] * apply(X.train, 2, sd) / sd(y.train))

# Testing Precision-Recall Curve
prc.df <- PrunePRC(vtm.test$threshold.df)

gg.prc <- ggplot(prc.df[nrow(prc.df):1, ], aes(recall, precision, color=threshold))
gg.prc <- ggPRC(gg.prc)
#prc.path <- file.path(graphics.dir, 'prc-testing.pdf')
#ggsave(prc.path, gg.prc, width=5, height=3.5)

# Testing ROC Curve
roc.df <- rbind(cbind(vtm.train$roc.df, 'part'='train'), cbind(vtm.test$roc.df, 'part'='test'))

gg.roc <- ggplot(roc.df, aes(fpr, recall, color=part))
gg.roc <- ggROC(gg.roc) + geom_path() +
  scale_color_manual(values=as.character(solarized[c('violet', 'green')]), 
    name='Partition (AUROC)', breaks=c('test', 'train'),
    labels=c(sprintf('Testing (%.3f)', vtm.test$auroc), 
             sprintf('Training (%.3f)', vtm.train$auroc)))
#roc.path <- file.path(graphics.dir, 'roc-testing.pdf')
#ggsave(roc.path, gg.roc, width=3, height=2.9)

# Save combined ROC and PRC
path <- file.path(graphics.dir, 'performance-testing.pdf')
OpenPDF(path, width=width.full, height=2.5)
gridExtra::grid.arrange(gg.roc, gg.prc, nrow=1, widths=c(1, 1.625))
ClosePDF(path)



