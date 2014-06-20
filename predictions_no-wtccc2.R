library(glmnet)
library(ggplot2)
library(doMC)
library(grid)
library(plyr)
library(reshape2)

options(stringsAsFactors=FALSE)

glmnet.alpha <- 0 # 0 for ridge regression

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
code.dir <- file.path(project.dir, 'rcode')
network.dir <- file.path(project.dir, 'networks', '140615-no-wtccc2')
feature.dir <- file.path(network.dir, 'features')

source(file.path(code.dir, 'machine-learning.R'))
source(file.path(code.dir, 'functions.R'))

feature.df <- ReadFeatures(feature.dir)
feature.names <- colnames(feature.df)[-(1:8)]

# Fit on whole data (global)
X <- as.matrix(feature.df[, feature.names])
y <- feature.df$status_int

X.model <- as.matrix(subset(feature.df, status_int != -1, select=feature.names))
y.model <- subset(feature.df, status_int != -1)[, 'status_int']


doMC::registerDoMC(cores=7)
cv.ridge <- glmnet::cv.glmnet(X.model, y.model, family='binomial', 
  alpha=glmnet.alpha, standardize=TRUE, parallel=TRUE)
# Save cv.glmnet object as an R data structure
model.path <- file.path(network.dir, 'model-ridge.rds')
saveRDS(cv.ridge, model.path)
cv.ridge <- readRDS(model.path)
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
coef.path <- file.path(network.dir, 'coefficients.txt')
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
  disease_name == 'multiple sclerosis' & status == 'HC_primary')$gene_symbol
genes.with.wtccc2 <- subset(associations.df.with.wtccc2, 
  disease_name == 'multiple sclerosis' & status == 'HC_primary')$gene_symbol
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


# Precision-Recall Curve
prc.df <- PrunePRC(vtm.wtc.novel$threshold.df)

gg.prc <- ggplot(prc.df[nrow(prc.df):1, ], aes(recall, precision, color=threshold))
gg.prc <- ggPRC(gg.prc) + scale_y_continuous(breaks=seq(0, 1, 0.1))

# ROC Curve
gg.roc <- ggplot(vtm.wtc.novel$roc.df, aes(fpr, recall))
gg.roc <- ggROC(gg.roc) + geom_path(color=as.character(solarized['violet'])) +
  annotate('text', x=Inf, y=-Inf, 
    label=sprintf('AUROC == %.3f', vtm.wtc.novel$auroc), 
    hjust=1.1, vjust=-1, parse=TRUE, size=4)

# Save combined ROC and PRC
perf.path <- file.path(network.dir, 'performance-wtccc2.pdf')
pdf(perf.path, width=width.full, height=2.5)
gridExtra::grid.arrange(gg.roc, gg.prc, nrow=1, widths=c(1, 1.625))
dev.off()





