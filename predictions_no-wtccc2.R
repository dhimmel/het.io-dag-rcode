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
network.dir <- file.path(project.dir, 'networks', '140522-no-wtccc2')
feature.dir <- file.path(network.dir, 'features')

source(file.path(code.dir, 'machine-learning.R'))
source(file.path(code.dir, 'functions.R'))

feature.df <- ReadFeatures(feature.dir)
feature.names <- colnames(feature.df)[-(1:8)]
# Exclude redundant MSigDB features
exclude.metanodes <- c('-C2.A-', '-C2.CP.A-', '-C3.A-', '-C4.A-', '-C5.A-')
exclude.feature <- rowSums(sapply(exclude.metanodes, grepl, feature.names, fixed=TRUE)) == 0
feature.names <- feature.names[exclude.feature]


# Fit on whole data (global)
X <- as.matrix(feature.df[, feature.names])
y <- feature.df$status_int

X.model <- as.matrix(subset(feature.df, status_int != -1, select=feature.names))
y.model <- subset(feature.df, status_int != -1)[, 'status_int']


doMC::registerDoMC(cores=7)
cv.ridge <- glmnet::cv.glmnet(X.model, y.model, family='binomial', 
  alpha=glmnet.alpha, standardize=TRUE, parallel=TRUE)
# Save cv.glmnet object as an R data structure
saveRDS(cv.ridge, file.path(network.dir, 'model.rds'))
lambda.global <- cv.ridge$lambda.1se

# Make predictions using the global model
## Save coefficients
coefs.global <- coef(cv.ridge, s=lambda.global)[, 1]
zcoefs.global <- c(0, coefs.global[-1] * apply(X.model, 2, sd) / sd(y.model))
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
feature.df[, 'prediction'] <- y.predicted
vtm <- VariableThresholdMetrics(
  subset(feature.df, status_int >= 0)$prediction, 
  subset(feature.df, status_int >= 0)$status_int)


# Read associated genes
associations.df.with.wtccc2 <- read.delim('/home/dhimmels/Documents/serg/data-sources/gwas-catalog/140205/processed/association-statuses.txt')
associations.df.no.wtccc2 <- read.delim('/home/dhimmels/Documents/serg/data-sources/gwas-catalog/140205/processed-no-wtccc2/association-statuses.txt')
genes.no.wtccc2 <- subset(associations.df.no.wtccc2, 
  disease_name == 'multiple sclerosis' & status == 'assoc_high')$gene_symbol
genes.with.wtccc2 <- subset(associations.df.with.wtccc2, 
  disease_name == 'multiple sclerosis' & status == 'assoc_high')$gene_symbol
genes.wtccc2.novel <- setdiff(genes.with.wtccc2, genes.no.wtccc2)
cat(sprintf('%s MS genes pre-WTCCC2. %s novel WTCCC2 genes.\n', length(genes.no.wtccc2), length(genes.wtccc2.novel)))

ms.df <- subset(feature.df, disease_name == 'multiple sclerosis')
ms.df$wtc_novel <- as.integer(ms.df$gene_symbol %in% genes.wtccc2.novel)
ms.features.path <- file.path(network.dir, 'MS-predictions-pre-wtccc2.txt')
write.table(ms.df, ms.features.path, sep='\t', row.names=FALSE, quote=FALSE)

ms.novel.df <- subset(ms.df, status_int != 1 | gene_symbol %in% genes.wtccc2.novel)
vtm.wtc.novel <- VariableThresholdMetrics(ms.novel.df$prediction, ms.novel.df$wtc_novel)

vtm.test.path <- file.path(network.dir, 'variable-threshold-metrics-testing.txt')
write.table(vtm.wtc.novel$threshold.df, vtm.test.path, sep='\t', row.names=FALSE, quote=FALSE)
#vtm.dapple <- VariableThresholdMetrics(ms.novel.df$dapple, ms.novel.df$wtc_novel)

# ROC Curve
roc.plot <- ggplot(vtm.wtc.novel$roc.df, aes(fpr, recall)) + 
  geom_abline(intercept=0, slope=1, color='red', linetype='dashed') +
  geom_line() + theme_bw() + coord_fixed() +
  #labels=c(sprintf('Testing (%.3f)', vtm.test$auc), sprintf('Training (%.3f)', vtm.train$auc))) +
  #theme(legend.justification=c(1,0), legend.position=c(1,0), legend.key=element_rect(linetype='blank')) +
  #geom_text(data=NULL, x=Inf, y=-Inf, label=sprintf('AUC=%.3f', vtm.hetnet$auc)) +
  annotate('text', x=Inf, y=-Inf, label=sprintf('AUC == %.3f', vtm.wtc.novel$auroc), hjust=1.1, vjust=-1, parse=TRUE) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) + 
  #theme(legend.background = element_rect(linetype='solid', color='grey')) +
  xlab('False Positive Rate') + ylab('Recall')
roc.path <- file.path(network.dir, 'roc.pdf')
ggsave(roc.path, roc.plot, width=3, height=2.9)


# Precision-Recall Curve
prc.plot <- ggplot(vtm.wtc.novel$threshold.df, aes(recall, precision, color=threshold)) + 
  geom_line(size=1, color='grey') + geom_point(size=1.5) + 
  theme_bw() + scale_color_gradientn(colours = rainbow(7)[1:6]) +
  theme(legend.justification=c(1, 1), legend.position=c(1, 1)) +
  theme(legend.background = element_rect(linetype='solid', color='grey')) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) +
  guides(color=guide_colorbar(title='Predicted\nprobability\nthreshold')) +
  xlab('Recall') + ylab('Precision')
prc.path <- file.path(network.dir, 'prc.pdf')
ggsave(prc.path, prc.plot, width=5, height=3.5)



