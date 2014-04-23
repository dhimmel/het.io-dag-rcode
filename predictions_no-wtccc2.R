library(glmnet)
library(ggplot2)
library(doMC)
library(grid)
library(plyr)
library(reshape2)

options(stringsAsFactors=FALSE)
options(width=Sys.getenv('COLUMNS'))

glmnet.alpha <- 0 # 0 for ridge regression

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
code.dir <- file.path(project.dir, 'rcode')
network.dir <- file.path(project.dir, 'networks', '140321-no-wtccc2')
feature.dir <- file.path(network.dir, 'features')

source(file.path(code.dir, 'machine-learning.R'))
feature.filenames <- list.files(feature.dir, pattern='DOID_')

feature.df.list <- list()
for (feature.filename in feature.filenames) {
  doid_code <- gsub('.txt.gz', '', feature.filename)
  cat(sprintf('Reading features for %s\n', doid_code))
  feature.path <- file.path(feature.dir, feature.filename)
  feature.df <- read.delim(feature.path, check.names=FALSE)
  # Remove genes that appear multiple times. Upstream HGNC bug
  feature.df <- feature.df[! duplicated(feature.df$source), ]
  feature.df.list[[doid_code]] <- feature.df
}

features.df <- do.call(rbind, feature.df.list)
feature.names <- colnames(feature.df)[-(1:5)]

doMC::registerDoMC(cores=7)
X <- as.matrix(features.df[, -(1:5)])
y <- features.df$status
cv.ridge <- glmnet::cv.glmnet(X, y, family='binomial', 
  alpha=glmnet.alpha, standardize=TRUE, parallel=TRUE)
# Save cv.glmnet object as an R data structure
saveRDS(cv.ridge, file.path(network.dir, 'model.rds'))
lambda.global <- cv.ridge$lambda.1se

# Make predictions using the global model
## Save coefficients
coefs.global <- coef(cv.ridge, s=lambda.global)[, 1]
zcoefs.global <- c(0, coefs.global[-1] * apply(X, 2, sd) / sd(y))
feature.name.mat <- do.call(rbind, strsplit(feature.names, '[|]'))
coef.df <- data.frame(
  'name'=c('intercept', feature.names), 
  'metric'=c(NA, feature.name.mat[, 1]), 
  'metapath'=c(NA, feature.name.mat[, 2]),
  'coef'=coefs.global,
  'zcoef'=zcoefs.global, row.names=NULL)
coef.path <- file.path(network.dir, 'coefs-global.txt')
write.table(coef.df, coef.path, sep='\t', row.names=FALSE, quote=FALSE)

## Make predictions
y.predicted <- predict(cv.ridge, s=lambda.global, newx=X, type='response')[, 1]
features.df[, 'prediction'] <- y.predicted
vtm.global <- VariableThresholdMetrics(features.df$prediction, features.df$status)


# Read associated genes
associations.df.with.wtccc2 <- read.delim('/home/dhimmels/Documents/serg/data-sources/gwas-catalog/140205/processed/associations.txt')
associations.df.no.wtccc2 <- read.delim('/home/dhimmels/Documents/serg/data-sources/gwas-catalog/140205/processed-no-wtccc2/associations.txt')
genes.no.wtccc2 <- subset(associations.df.no.wtccc2, doid_name == 'multiple sclerosis')$symbol
genes.with.wtccc2 <- subset(associations.df.with.wtccc2, doid_name == 'multiple sclerosis')$symbol
genes.wtccc2.novel <- setdiff(genes.with.wtccc2, genes.no.wtccc2)

cat(sprintf('%s MS genes pre-WTCCC2. %s novel WTCCC2 genes.\n', length(genes.no.wtccc2), length(genes.wtccc2.novel)))

ms.df <- subset(features.df, target_name == 'multiple sclerosis')
ms.df$wtc_novel <- as.integer(ms.df$source %in% genes.wtccc2.novel)
ms.features.path <- file.path(network.dir, 'MS-predictions-pre-wtccc2.txt')
write.table(ms.df, ms.features.path, sep='\t', row.names=FALSE, quote=FALSE)

ms.novel.df <- subset(ms.df, status == 0)
vtm.hetnet <- VariableThresholdMetrics(ms.novel.df$prediction, ms.novel.df$wtc_novel)
vtm.test.path <- file.path(network.dir, 'variable-threshold-metrics-testing.txt')
write.table(vtm.hetnet$threshold.df, vtm.test.path, sep='\t', row.names=FALSE, quote=FALSE)
#vtm.dapple <- VariableThresholdMetrics(ms.novel.df$dapple, ms.novel.df$wtc_novel)

# ROC Curve
roc.plot <- ggplot(vtm.hetnet$threshold.df, aes(fpr, tpr)) + 
  geom_abline(intercept=0, slope=1, color='red', linetype='dashed') +
  geom_line() + theme_bw() + coord_fixed() +
  #labels=c(sprintf('Testing (%.3f)', vtm.test$auc), sprintf('Training (%.3f)', vtm.train$auc))) +
  #theme(legend.justification=c(1,0), legend.position=c(1,0), legend.key=element_rect(linetype='blank')) +
  #geom_text(data=NULL, x=Inf, y=-Inf, label=sprintf('AUC=%.3f', vtm.hetnet$auc)) +
  annotate('text', x=Inf, y=-Inf, label=sprintf('AUC == %.3f', vtm.hetnet$auc), hjust=1.1, vjust=-1, parse=TRUE) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) + 
  #theme(legend.background = element_rect(linetype='solid', color='grey')) +
  xlab('False Positive Rate') + ylab('Recall')
roc.path <- file.path(network.dir, 'roc.pdf')
ggsave(roc.path, roc.plot, width=3, height=2.9)


# Precision-Recall Curve
prc.plot <- ggplot(vtm.hetnet$threshold.df, aes(recall, precision, color=threshold)) + 
  geom_line(size=1, color='grey') + geom_point(size=1.5) + 
  theme_bw() + scale_color_gradientn(colours = rainbow(7)[1:6]) +
  theme(legend.justification=c(1, 1), legend.position=c(1, 1)) +
  theme(legend.background = element_rect(linetype='solid', color='grey')) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) +
  guides(color=guide_colorbar(title='Predicted\nprobability\nthreshold')) +
  xlab('Recall') + ylab('Precision')
prc.path <- file.path(network.dir, 'prc.pdf')
ggsave(prc.path, prc.plot, width=5, height=3.5)



