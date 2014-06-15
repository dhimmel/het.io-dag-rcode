library(glmnet)
library(ggplot2)
library(doMC)
library(grid)
library(gtools)
library(plyr)
library(reshape2)

options(stringsAsFactors=FALSE)
options(width=Sys.getenv('COLUMNS'))

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
network.id <- '140614-training'

code.dir <- file.path(project.dir, 'rcode')
network.dir <- file.path(project.dir, 'networks', network.id)
feature.dir <- file.path(network.dir, 'features')
graphics.dir <- file.path(network.dir, 'graphics')
modeling.dir <- file.path(network.dir, 'modeling')
dir.create(graphics.dir, showWarnings=FALSE)
dir.create(modeling.dir, showWarnings=FALSE)

source(file.path(code.dir, 'machine-learning.R'))
source(file.path(code.dir, 'functions.R'))

feature.df <- ReadFeatures(feature.dir)
feature.names <- colnames(features.df)[-(1:8)]



glmnet.alpha <- 0 # 0 for ridge regression

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
code.dir <- file.path(project.dir, 'rcode')
network.dir <- file.path(project.dir, 'networks', '140522-training')
feature.dir <- file.path(network.dir, 'features')
graphics.dir <- file.path(network.dir, 'graphics')
modeling.dir <- file.path(network.dir, 'modeling')
webdata.dir <- file.path(network.dir, 'webdata')
dir.create(graphics.dir, showWarnings=FALSE)
dir.create(modeling.dir, showWarnings=FALSE)
dir.create(webdata.dir, showWarnings=FALSE)

# Read Features
feature.df <- ReadFeatures(feature.dir)
feature.names <- colnames(feature.df)[-(1:8)]

train.df <- subset(feature.df, part=='train' & status_int != -1)
test.df <- subset(feature.df, part=='test' & status_int != -1)
#train.df <- subset(feature.df, percentile > 0.20 & status_int != -1)
#test.df <- subset(feature.df, percentile <= 0.20 & status_int != -1)

X.train <- as.matrix(train.df[, feature.names])
X.test <- as.matrix(test.df[, feature.names])
y.train <- train.df$status_int
y.test <- test.df$status_int

doMC::registerDoMC(cores=7)
set.seed(0)
cv.ridge.train <- glmnet::cv.glmnet(X.train, y.train, family='binomial', 
  alpha=glmnet.alpha, standardize=TRUE, parallel=TRUE)
lambda.train <- cv.ridge.train$lambda.1se
# Save cv.glmnet object as an R data structure
saveRDS(cv.ridge.train, file.path(modeling.dir, 'training-model.rds'))
#cv.ridge.train <- readRDS(file.path(network.dir, 'training-model.rds'))
y.predicted.train <- as.numeric(predict(cv.ridge.train, s=lambda.train, newx=X.train, type='response'))
y.predicted.test <- as.numeric(predict(cv.ridge.train, s=lambda.train, newx=X.test, type='response'))
vtm.train <- VariableThresholdMetrics(y.predicted.train, y.train)
vtm.test <- VariableThresholdMetrics(y.predicted.test, y.test)

vtm.train.df <- vtm.train$threshold.df
vtm.test.df <- vtm.test$threshold.df
vtm.part.df <- rbind(cbind(vtm.train.df, 'part'='train'), cbind(vtm.test.df, 'part'='test'))

vtm.test.path <- file.path(modeling.dir, 'variable-threshold-metrics-testing.txt')
write.table(round(vtm.test.df, 5), vtm.test.path, sep='\t', row.names=FALSE, quote=FALSE)


coefs.training <- coef(cv.ridge.train, s=lambda.train)[, 1]
zcoefs.training <- c(0, coefs.training[-1] * apply(X.train, 2, sd) / sd(y.train))

# Testing Precision-Recall Curve
prc.plot <- ggplot(vtm.test.df[nrow(vtm.test.df):1, ], aes(recall, precision, color=threshold)) + 
  geom_line(size=1, color='grey') + geom_point(size=1.5) + 
  theme_bw() + scale_color_gradientn(colours = rainbow(7)[1:6]) +
  theme(legend.justification=c(1, 1), legend.position=c(1, 1)) +
  theme(legend.background = element_rect(linetype='solid', color='grey')) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) +
  guides(color=guide_colorbar(title='Predicted\nprobability\nthreshold')) +
  xlab('Recall') + ylab('Precision')
prc.path <- file.path(graphics.dir, 'prc-testing.pdf')
ggsave(prc.path, prc.plot, width=5, height=3.5)

# Testing ROC Curve
roc.df <- rbind(cbind(vtm.train$roc.df, 'part'='train'), cbind(vtm.test$roc.df, 'part'='test'))
roc.plot <- ggplot(roc.df, aes(fpr, recall, linetype=part)) + 
  geom_line() + theme_bw() + coord_fixed() +
  scale_linetype_manual(values=c('solid', 'dashed'), 
    name='Partition (AUC)', breaks=c('test', 'train'),
    labels=c(sprintf('Testing (%.3f)', vtm.test$auroc), sprintf('Training (%.3f)', vtm.train$auroc))) +
  theme(legend.justification=c(1,0), legend.position=c(1,0), legend.key=element_rect(linetype='blank')) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) + 
  theme(legend.background = element_rect(linetype='solid', color='grey')) +
  xlab('False Positive Rate') + ylab('Recall')
roc.path <- file.path(graphics.dir, 'roc-testing.pdf')
ggsave(roc.path, roc.plot, width=3, height=2.9)


