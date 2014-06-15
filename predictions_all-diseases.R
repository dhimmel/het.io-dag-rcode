library(glmnet)
library(ggplot2)
library(doMC)
library(grid)
library(gtools)
library(plyr)
library(reshape2)

options(stringsAsFactors=FALSE)
options(width=Sys.getenv('COLUMNS'))

glmnet.alpha <- 0 # 0 for ridge regression


project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
network.id <- '140522-all-assoc-lessmsig'

code.dir <- file.path(project.dir, 'rcode')
network.dir <- file.path(project.dir, 'networks', network.id)
feature.dir <- file.path(network.dir, 'features')
graphics.dir <- file.path(network.dir, 'graphics')
modeling.dir <- file.path(network.dir, 'modeling')
webdata.dir <- file.path(network.dir, 'webdata')
dir.create(graphics.dir, showWarnings=FALSE)
dir.create(modeling.dir, showWarnings=FALSE)
dir.create(webdata.dir, showWarnings=FALSE)

source(file.path(code.dir, 'machine-learning.R'))
source(file.path(code.dir, 'functions.R'))

features.df <- ReadFeatures(feature.dir)
feature.names <- colnames(features.df)[-(1:8)]

# Exclude redundant MSigDB features
exclude.metanodes <- c('-C2.A-', '-C2.CP.A-', '-C3.A-', '-C4.A-', '-C5.A-')
exclude.feature <- rowSums(sapply(exclude.metanodes, grepl, feature.names, fixed=TRUE)) == 0
feature.names <- feature.names[exclude.feature]


category.path <- file.path(project.dir, 'data-integration', 'disease-categories.txt')
category.df <- read.delim(category.path)
features.df$disease_category <- category.df[match(features.df$disease_code, category.df$doid_code), 'category']

pathophys.path <- file.path(project.dir, 'data-integration', 'pathophysiology.txt')
pathophys.df <- read.delim(pathophys.path)
features.df$disease_pathophys <- pathophys.df[match(features.df$disease_code, pathophys.df$disease_code), 'pathophysiology']

################################################################################
## Fit on all oberservations

# Fit on whole data (global)
X <- as.matrix(features.df[, feature.names])
y <- features.df$status_int

X.model <- as.matrix(subset(features.df, status_int != -1, select=feature.names))
y.model <- subset(features.df, status_int != -1)[, 'status_int']

doMC::registerDoMC(cores=7)
set.seed(0)
cv.ridge <- glmnet::cv.glmnet(X.model, y.model, family='binomial', 
  alpha=glmnet.alpha, standardize=TRUE, parallel=TRUE)
# Save cv.glmnet object as an R data structure
saveRDS(cv.ridge, file.path(modeling.dir, 'global-model.rds'))
#cv.ridge <- readRDS(file.path(modeling.dir, 'global-model.rds'))
lambda.global <- cv.ridge$lambda.1se

# Make predictions using the global model
## Save coefficients
coefs.global <- coef(cv.ridge, s=lambda.global)[, 1]
zcoefs.global <- c(0, coefs.global[-1] * apply(X.model, 2, sd) / sd(y.model))
feature.name.mat <- do.call(rbind, strsplit(feature.names, '[|]'))
coef.df <- data.frame(
  'feature'=c('intercept', feature.names), 
  'metric'=c(NA, feature.name.mat[, 1]), 
  'metapath'=c(NA, feature.name.mat[, 2]),
  'ridge_coef'=coefs.global,
  'ridge_zcoef'=zcoefs.global, row.names=NULL)

## Save XB and make predictions
XB.mat <- X %*% diag(coefs.global[-1])
XB.names <- paste('XB|', feature.names, sep='')
colnames(XB.mat) <- XB.names
y.predicted <- as.numeric(predict(cv.ridge, s=lambda.global, newx=X, type='response'))
features.df[, 'prediction'] <- y.predicted
prediction.df <- cbind(features.df, as.data.frame(XB.mat, check.names=FALSE))
predictions.file <- gzfile(file.path(modeling.dir, 'predictions.txt.gz'), 'w')
write.table(prediction.df, predictions.file, sep='\t', row.names=FALSE, quote=FALSE)
close(predictions.file)
feat.perf.df <- subset(features.df, status_int != -1)
vtm.global <- VariableThresholdMetrics(feat.perf.df$prediction, feat.perf.df$status_int)

prediction.cast <- reshape2::dcast(features.df, gene_symbol ~ disease_name, value.var='prediction')
predictions.cast.path <- file.path(modeling.dir, 'prediction-table.txt')
write.table(prediction.cast, predictions.cast.path, sep='\t', row.names=FALSE, quote=FALSE)


################################################################################
## Parsimonious Models

set.seed(0)
cv.lasso <- glmnet::cv.glmnet(X.model, y.model, family='binomial', 
  alpha=1, standardize=TRUE, parallel=TRUE)

lasso.predicted <- as.numeric(predict(cv.ridge, s=cv.lasso$lambda.1se, newx=X.model, type='response'))
lasso.vtm <- VariableThresholdMetrics(lasso.predicted, feat.perf.df$status_int)


coef.df$lasso_coef <- coef(cv.lasso, s=cv.lasso$lambda.1se)[, 1]
coef.df$lasso_zcoef <- c(0, coef.df$lasso_coef[-1] * apply(X.model, 2, sd) / sd(y.model))

coef.path <- file.path(modeling.dir, 'coefs-global.txt')
write.table(coef.df, coef.path, sep='\t', row.names=FALSE, quote=FALSE)


feature.desc.path <- file.path(project.dir, 'data-integration', 'feature-descriptions.txt')
desc.df <- read.delim(feature.desc.path)

feature.converter <- desc.df$easy_name
names(feature.converter) <- desc.df$feature

coef.df$neg_ridge_zcoef <- -coef.df$ridge_zcoef
coef.melt <- reshape2::melt(subset(coef.df, feature != 'intercept'),
  measure.vars=c('neg_ridge_zcoef', 'lasso_zcoef'))

coef.melt.ridge <- subset(coef.melt, variable=='neg_ridge_zcoef')
features.sorted <- coef.melt.ridge[order(coef.melt.ridge$value, decreasing=TRUE), 'feature']


fill.red <- '#FB7C72'
fill.blue <- '#7284FF'

gg.zcoef <- ggplot(coef.melt, aes(x=feature, ymin=0, ymax=value, color=variable)) + 
  theme_bw() + geom_hline(yintercept=0) + ylab('Standardized Coefficient') +
  geom_linerange(stat='identity', size=5) + coord_flip() +
  scale_x_discrete(limits=features.sorted, labels=feature.converter[features.sorted]) + 
  scale_y_continuous(labels=abs, limits=c(-max(coef.melt$value), max(coef.melt$value))) +
  scale_color_manual(name='Method (AUROC)', values=c(fill.red, fill.blue), 
    labels=c(sprintf('ridge (%.3f)', vtm.global$auroc), sprintf('lasso (%.3f)', lasso.vtm$auroc))) +
  theme(legend.key=element_rect(linetype='blank')) +
  theme(axis.text.y=element_text(angle=30, hjust=1)) + xlab(NULL) +
  theme(legend.justification=c(1,0), legend.position=c(1,0), 
    legend.background=element_rect(color='grey60', size=0.25))
ggsave(file.path(graphics.dir, 'coefficients.pdf'), gg.zcoef, width=4.5, height=5)


################################################################################
## Performance by status type
violin.plot <- ggplot(features.df, aes(status, prediction)) +
  geom_violin() + coord_cartesian(ylim=c(0.0008, 0.01)) + 
  scale_y_log10() + stat_summary(fun.y=mean, geom='point', color='darkgreen', size=3) +
  theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(file.path(graphics.dir, 'prediction-violins.pdf'), violin.plot, width=4.5, height=5)

pos.statuses <- c('assoc_high', 'linked_high', 'assoc_low', 'linked_low')
auroc.list <- list()

roc.df <- do.call(rbind, lapply(pos.statuses, function(pos.status) {
  subset.df <- subset(features.df, status %in% c('negative', pos.status))
  status <- subset.df$status == pos.status
  vtm <- VariableThresholdMetrics(subset.df$prediction, status)
  roc.df <- vtm$roc.df
  roc.df[, 'Positives'] <- pos.status
  #split.status <- strsplit(pos.status, '_')[[1]]
  #roc.df[, 'Association'] <- split.status[1]
  #roc.df[, 'Confidence'] <- split.status[2]
  auroc.list[[pos.status]] <<- vtm$auroc
  return(roc.df)
}))


cols <- as.character(solarized[c('red', 'blue', 'red', 'blue')])
linetypes <- c('solid', 'solid', 'dotted', 'dotted')
gglabels <- c(
  sprintf('HC Primary (%.2f)', auroc.list$assoc_high),
  sprintf('HC Linked (%.2f)', auroc.list$linked_high),
  sprintf('LC Primary (%.2f)', auroc.list$assoc_low),
  sprintf('LC Linked (%.2f)', auroc.list$linked_low))

gg.roc <- ggplot(roc.df, aes(fpr, recall, color=Positives, linetype=Positives))
gg.roc <- ggROC(gg.roc) + geom_line(size=0.9) +
  scale_linetype_manual(name='Positive Set (AUROC)', values=linetypes, labels=gglabels, breaks=pos.statuses) +
  scale_color_manual(name='Positive Set (AUROC)', values=cols, labels=gglabels, breaks=pos.statuses)
ggsave(file.path(graphics.dir, 'ROC-by-positive-set.pdf'), gg.roc, width=3.7, height=3.7)


################################################################################
## Feature correlation plot
HclustOrder <- function(feature.mat) {
  # Returns the colnames, ordered by heirarchical clustering
  cor.mat <- cor(feature.mat)
  clust <- hclust(dist(cor.mat), method='ward')
  return(colnames(feature.mat)[clust$order])
}

features.sorted <- HclustOrder(X)
cor.mat <- X[, features.sorted]
cor.tri.mat <- cor(cor.mat)
cor.tri.mat[lower.tri(cor.tri.mat, diag=TRUE)] <- NA
cor.melt <- na.omit(as.data.frame(reshape2::melt(cor.tri.mat, varnames=c('x', 'y'), value.name='correlation')))


xlimits <- features.sorted[-length(features.sorted)]
ylimits <- rev(features.sorted[-1])

diverging.cols <- c('#0000FF', '#f7eff7', '#FF0000') # blues
cor.plot <- ggplot(cor.melt, aes(x, y, fill=correlation)) +
  geom_tile(color='white') + 
  scale_fill_gradientn(name=expression("Pearson's" * ~ rho), limits=c(-1, 1), colours=diverging.cols) + 
  #scale_fill_gradient2(name=expression("Pearson's" * ~ rho), limits=c(-1, 1)) + 
  scale_x_discrete(limits=xlimits, labels=feature.converter[xlimits]) + 
  scale_y_discrete(limits=ylimits, labels=feature.converter[ylimits]) +
  theme_bw() + theme(axis.text.x=element_text(angle=65, hjust=1)) +
  xlab(NULL) + ylab(NULL) + coord_fixed() +
  theme(plot.background=element_blank(), panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), panel.border=element_blank(),
    axis.line=element_line(color='black')) +
  theme(legend.justification=c(1, 0), 
    legend.position=c(1.0, 0.75), legend.direction='horizontal') +
  guides(fill=guide_colorbar(barwidth=10, barheight=1, title.position='top', title.hjust=0.5))

ggsave(file.path(graphics.dir, 'feature-correlations.pdf'), cor.plot, width=7, height=7)


################################################################################
## Disease Specific Predictions on global data set

# auroc.df -- data.frame of AUC under different conditions
# breadth (disease_specific or global), type (feature or model), feature, metric, metapath, name, disease_code, disease_name, disease_category, positives, auroc
auroc.colnames <- c('breadth', 'type', 'name', 'feature', 'metric', 'metapath', 'disease_code', 'disease_name', 'disease_category', 'positives', 'auroc')
auroc.blank.df <- data.frame(t(rep(NA, length(auroc.colnames))))
names(auroc.blank.df) <- auroc.colnames
auroc.blank.df <- auroc.blank.df[-1, ]

OrderAurocDf <- function(auroc.df) {
  auroc.df <- merge(auroc.df, auroc.blank.df, all.x=TRUE, all.y=TRUE)
  auroc.df <- auroc.df[, match(auroc.colnames, colnames(auroc.df))]
  return(auroc.df)
}


# Calculate disease specific AUCs from the global ridge model
auroc.df1 <- plyr::ddply(feat.perf.df, c('disease_code'), plyr::summarize,
  'breadth' = 'disease_specific',
  'type'='model', 'name'='ridge',
  'auroc' = VariableThresholdMetrics(prediction, status_int)$auroc,
  'positives'=sum(status_int))
auroc.df1 <- OrderAurocDf(auroc.df1)

ComputeFeatureAUC <- function(feat.df) {
  apply(feat.df[, feature.names], 2, function(x) VariableThresholdMetrics(x, feat.df$status_int)$auroc)
}

# Calculate global AUCs for each feature
auroc.vec2 <- ComputeFeatureAUC(feat.perf.df)
auroc.df2 <- data.frame(
  'feature'=names(auroc.vec2), 
  'auroc'=as.numeric(auroc.vec2),
  'type'='feature',
  'breadth'='global',
  'positives'=sum(feat.perf.df$status_int)
)
auroc.df2 <- OrderAurocDf(auroc.df2)

# Calculate disease-specific AUROCs for each feature
fXd.auroc.tab <- plyr::ddply(feat.perf.df, c('disease_code'), ComputeFeatureAUC)
auroc.df3 <- reshape2::melt(fXd.auroc.tab, id.vars=c('disease_code'), variable.name='feature', value.name='auroc')
auroc.df3 <- cbind(auroc.df3, 'type'='feature', 'breadth'='disease_specific')
auroc.df3$positives <- category.df[match(auroc.df3$disease_code, category.df$doid_code), 'count']
auroc.df3 <- OrderAurocDf(auroc.df3)

# Global AUC for Ridge Model
auroc.df4 <- data.frame('breadth'='global', 'type'='model', 'name'='ridge',
  'positives'=sum(feat.perf.df$status_int), 'auroc'=vtm.global$auroc)
auroc.df4 <- OrderAurocDf(auroc.df4)

# Combine
auroc.df <- rbind(auroc.df1, auroc.df2, auroc.df3, auroc.df4)

auroc.df$disease_category <- category.df[match(auroc.df$disease_code, category.df$doid_code), 'category']
auroc.df$disease_pathophys <- pathophys.df[match(auroc.df$disease_code, pathophys.df$disease_code), 'pathophysiology']
auroc.df$disease_name <- pathophys.df[match(auroc.df$disease_code, pathophys.df$disease_code), 'disease_name']

metric.metapath.mat <- do.call(rbind, strsplit(auroc.df$feature, '|', fixed=TRUE))
auroc.df$metric <- metric.metapath.mat[, 1]
auroc.df$metapath <- gsub('-', '', metric.metapath.mat[, 2])
is.feature.auc <- ! is.na(auroc.df$feature)
auroc.df[is.feature.auc, 'name'] <- auroc.df[is.feature.auc, 'metapath']


################################################################################
## Plot AUROCs

# Read MSigDB Nomenclature file
msig.nomen.path <- file.path(project.dir, 'data-integration', 'MSigDB-nomenclature.txt')
msig.nomen.df <- read.delim(msig.nomen.path)
msig.nomen.df$metapath.abbrev <- paste('Gm', msig.nomen.df$code, 'mGaD', sep='')

NAtoFALSE <- function(vec) {
  vec[is.na(vec)] <- FALSE
  return(vec)
}
auroc.df$panel <- 'Model'
auroc.df[NAtoFALSE(substr(auroc.df$metric, 1, 4)) == 'PC_s', 'panel'] <- 'PCs'
auroc.df[NAtoFALSE(substr(auroc.df$metric, 1, 4)) == 'PC_t', 'panel'] <- 'PCt'
auroc.df[NAtoFALSE(substr(auroc.df$metric, 1, 4)) == 'DWPC', 'panel'] <- 'DWPC'
is.msig.auc <- auroc.df$metapath %in% msig.nomen.df$metapath.abbrev
auroc.df[is.msig.auc, 'panel'] <- 'Gene-{MSigDB Collection}-Gene-Disease DWPC Feature'
auroc.df[is.msig.auc, 'name'] <- msig.nomen.df[match(auroc.df[is.msig.auc, 'metapath'], msig.nomen.df$metapath.abbrev), 'short_name']

MeanConfInt <- function(x) {t.test(x)$conf.int[1:2]}
category.colors <- c('#005200', '#B20000', '#8F008F', '#0000B2', '#E68A00')
pathophys.colors <- c('#005200', '#B20000', '#8F008F', '#0000B2', '#E68A00', 'black')

#pathophys.colors <- as.character(solarized[c('red', 'yellow', 'magenta', 'violet', 'cyan', 'base01')])

## MSigDB Plot
msig.df <- subset(auroc.df, panel == 'Gene-{MSigDB Collection}-Gene-Disease DWPC Feature' | panel == 'Model')
msig.disease.df <- subset(msig.df, breadth == 'disease_specific')
msig.global.df <- subset(msig.df, breadth == 'global')
# Order by global AUROC
msig.levels <- as.character(msig.global.df$name[order(msig.global.df$auroc)])
msig.disease.df$name <- factor(msig.disease.df$name, levels=msig.levels)
msig.global.df$name <- factor(msig.global.df$name, levels=msig.levels)

set.seed(0); msig.plot <- ggplot(msig.disease.df, aes(name, auroc)) + 
  facet_grid(. ~ panel, scales='free_x', space='free_x') +
  geom_hline(aes(yintercept=0.5), color='darkgrey', linetype='dashed') +
  stat_summary(fun.y='MeanConfInt', geom='line', color='grey', size=9.5) +
  geom_point(data=msig.global.df, size=11.5, shape='-', color='#4C4C4C') + 
  geom_point(aes(color=disease_pathophys), position=position_jitter(width=0.28), alpha=0.7, size=2.5) + 
  theme(plot.margin = grid::unit(c(0, 0, 0, 0), 'points')) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  xlab(NULL) + ylab('AUROC') +# guides(color=FALSE) +
  theme(legend.position='top', legend.margin=grid::unit(-10, 'points')) +
  #theme(legend.justification=c(0,1), legend.position=c(0, 1),
  #  legend.background=element_rect(color='grey60', size=0.25),
  #  legend.margin=grid::unit(1, 'points'), legend.direction='horizontal') +
  scale_colour_manual(name='Pathophysiology', values=pathophys.colors) +
  guides(color=guide_legend(nrow=2))
msig.plot.path <- file.path(graphics.dir, 'AUROC-msig-plot.pdf')
ggsave(msig.plot.path, msig.plot, width=6.83, height=5)

## NonMSigDB Plot
nonmsig.df <- subset(auroc.df, panel != 'Gene-{MSigDB Collection}-Gene-Disease DWPC Feature')
nonmsig.df$panel <- factor(nonmsig.df$panel, levels=c('DWPC', 'PCt', 'PCs', 'Model'))
nonmsig.disease.df <- subset(nonmsig.df, breadth == 'disease_specific' & panel != 'PCt')
nonmsig.global.df <- subset(nonmsig.df, breadth == 'global')

# Order by global AUROC
nonmsig.levels <- unique(as.character(nonmsig.global.df$name[order(nonmsig.global.df$auroc)]))
nonmsig.disease.df$name <- factor(nonmsig.disease.df$name, levels=nonmsig.levels)
nonmsig.global.df$name <- factor(nonmsig.global.df$name, levels=nonmsig.levels)

set.seed(0); nonmsig.plot <- ggplot(nonmsig.disease.df, aes(name, auroc)) + 
  facet_grid(. ~ panel, scales='free_x', space='free_x') +
  geom_hline(aes(yintercept=0.5), color='darkgrey', linetype='dashed') +
  stat_summary(fun.y='MeanConfInt', geom='line', color='grey', size=12) +
  geom_point(data=nonmsig.global.df, size=15, shape='-', color='#4C4C4C') + 
  geom_point(aes(color=disease_pathophys), position=position_jitter(width=0.3), alpha=0.7, size=2.5) + 
  theme(plot.margin = grid::unit(c(0, 0, 0, 0), 'points')) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  xlab(NULL) + ylab('AUROC') + 
  theme(legend.key=element_rect(linetype='blank')) +
  #theme(legend.position='top', legend.margin=grid::unit(-10, 'points')) +
  #theme(legend.justification=c(0,1), legend.position=c(0, 1),
  #  legend.background=element_rect(color='grey60', size=0.25),
  #  legend.margin=grid::unit(1, 'points'), legend.direction='horizontal') +
  theme(legend.margin=grid::unit(1, 'points')) +
  scale_colour_manual(name='Pathophysiology', values=pathophys.colors)
nonmsig.plot.path <- file.path(graphics.dir, 'AUROC-nomsig-plot.pdf')
ggsave(nonmsig.plot.path, nonmsig.plot, width=6.83, height=4)


################################################################################
# Webdata

#hgnc.path <- '/home/dhimmels/Documents/serg/data-sources/hgnc/140205/protein-coding.txt'
#hgnc.df <- read.delim(hgnc.path)
#features.df$gene_code <- hgnc.df[match(features.df$gene_symbol, hgnc.df$symbol), 'hgnc_id']

#################
# Disease Summary Table
# disease_name, disease_code, disease_category, associations, auroc
disease.summary.df <- subset(auroc.df, breadth == 'disease_specific' & type == 'model',
  select=c('disease_name', 'disease_code', 'disease_pathophys', 'positives', 'auroc'))

colnames(disease.summary.df)[colnames(disease.summary.df) == 'positives'] <- 'associations'
colnames(disease.summary.df)[colnames(disease.summary.df) == 'disease_pathophys'] <- 'pathophysiology'
disease.summary.path <- file.path(webdata.dir, 'disease-summary-table.txt')
disease.summary.df <- format(disease.summary.df, digits=2)
write.table(disease.summary.df, disease.summary.path, sep='\t', row.names=FALSE, quote=FALSE)

#################
# Gene Summary Table
#gene_symbol, gene_code, associations, mean_prediction
gene.summary.df <- plyr::ddply(features.df, c('gene_symbol', 'gene_code'), plyr::summarize,
  'associations'=sum(status_int == 1),
  'mean_prediction'=mean(prediction)
)
gene.summary.path <- file.path(webdata.dir, 'gene-summary-table.txt')
gene.summary.df <- format(gene.summary.df, digits=2)
write.table(gene.summary.df, gene.summary.path, sep='\t', row.names=FALSE, quote=FALSE)

#################
# Feature Summary Table

# metric, metapath, zcoef, coef, auroc;
feature.summary.df <- subset(auroc.df, breadth == 'global' & type == 'feature',
  select=c('feature', 'name', 'metric', 'auroc'))
colnames(feature.summary.df)[colnames(feature.summary.df) == 'name'] <- 'metapath'
feature.summary.df[feature.summary.df$metric == 'PC_s', 'metric'] <- 'PC_source'
feature.summary.df[feature.summary.df$metric == 'PC_t', 'metric'] <- 'PC_target'
feature.summary.df$standardized_coefficient <- coef.df[match(feature.summary.df$feature, coef.df$feature), 'ridge_zcoef']
feature.summary.df$feature <- gsub('|', '_', feature.summary.df$feature, fixed=TRUE)
feature.summary.path <- file.path(webdata.dir, 'feature-summary-table.txt')
feature.summary.df <- format(feature.summary.df, digits=2)
write.table(feature.summary.df, feature.summary.path, sep='\t', row.names=FALSE, quote=FALSE)


#################
# Disease Tables

status.converter <- c('assoc_high'='+ HCA', 'assoc_low'='\u00B1 LCA', 'linked_high'='\u00B1 HCL', 'linked_low'='\u00B1 LCL', 'negative'='-')

dir.create(file.path(webdata.dir, 'disease-tables'), showWarnings=FALSE)
#gene_symbol, gene_code, positives, mean_prediction; prediction
MakeDiseaseTable <- function(feature.df) {
  gene.summary.index <- match(feature.df$gene_symbol, gene.summary.df$gene_symbol)
  disease.table <- data.frame(
    'gene_symbol'=feature.df$gene_symbol, 
    'gene_code'=feature.df$gene_code, 
    'status'=status.converter[feature.df$status], 
    'other_associations'=feature.df[, 'PC_s|G-a-D'],
    'mean_prediction'=gene.summary.df[gene.summary.index, 'mean_prediction'],
    'prediction'=feature.df$prediction)
  return(disease.table)
}
disease.tables <- plyr::dlply(features.df, 'disease_code', MakeDiseaseTable)
for (disease_code in names(disease.tables)) {
  disease.table <- disease.tables[[disease_code]]
  disease.table <- format(disease.table, digits=2)
  filename <- sprintf('%s.txt', gsub(':', '_', disease_code, fixed=TRUE))
  table.path <- file.path(webdata.dir, 'disease-tables', filename)
  write.table(disease.table, table.path, sep='\t', row.names=FALSE, quote=FALSE)
}

#################
# Gene Tables
dir.create(file.path(webdata.dir, 'gene-tables'), showWarnings=FALSE)
#disease_name; doid_code; category; associations; mean_prediction; prediction
MakeGeneTable <- function(feature.df) {
  #disease.summary.index <- match(feature.df$disease_code, disease.summary.df$disease_code)
  gene.table <- data.frame(
    'disease_name'=feature.df$disease_name, 
    'disease_code'=feature.df$disease_code, 
    'pathophysiology'=feature.df$disease_pathophys,
    'status'=status.converter[feature.df$status], 
    'other_associations'=feature.df[, 'PC_t|G-a-D'],
    'prediction'=feature.df$prediction)
  return(gene.table)
}
gene.tables <- plyr::dlply(features.df, 'gene_code', MakeGeneTable)
for (gene.code in names(gene.tables)) {
  gene.table <- gene.tables[[gene.code]]
  gene.table <- format(gene.table, digits=2)
  filename <- sprintf('%s.txt', gsub(':', '_', gene.code, fixed=TRUE))
  table.path <- file.path(webdata.dir, 'gene-tables', filename)
  write.table(gene.table, table.path, sep='\t', row.names=FALSE, quote=FALSE)
}


#################
# Feature Tables
dir.create(file.path(webdata.dir, 'feature-tables'), showWarnings=FALSE)
#disease_name, doid_code; associations; global_auroc, feature_auroc
fXd.auroc.df <- subset(auroc.df, breadth == 'disease_specific' & type == 'feature')
fXd.auroc.df$feature <- gsub('|', '_', fXd.auroc.df$feature, fixed=TRUE)

MakeFeatureTable <- function(feature.auroc.df) {
  feature.table <- data.frame(
    'disease_name'=feature.auroc.df$disease_name,
    'disease_code'=feature.auroc.df$disease_code,
    'pathophysiology'=feature.auroc.df$disease_pathophys,
    'associations'=feature.auroc.df$positives,
    'disease_global_auroc'=disease.summary.df[match(feature.auroc.df$disease_code, disease.summary.df$disease_code), 'auroc'],
    'auroc'=feature.auroc.df$auroc)
  return(feature.table)
}

feature.tables <- plyr::dlply(fXd.auroc.df, 'feature', MakeFeatureTable)
for (feature in names(feature.tables)) {
  feature.table <- feature.tables[[feature]]
  feature.table <- format(feature.table, digits=2)
  filename <- sprintf('%s.txt', feature)
  table.path <- file.path(webdata.dir, 'feature-tables', filename)
  write.table(feature.table, table.path, sep='\t', row.names=FALSE, quote=FALSE)
}


#################
# Coef Tables


coef.match.name.df <- subset(auroc.df, breadth == 'global' & type == 'feature',
  select=c('feature', 'name', 'metric', 'panel'))
coef.match.name.df$easy_name <- coef.match.name.df$name
coef.match.name.df[coef.match.name.df$feature == 'PC_s|G-a-D', 'easy_name'] <- 'GaD|PC_source'
coef.match.name.df[coef.match.name.df$feature == 'PC_t|G-a-D', 'easy_name'] <- 'GaD|PC_target'
msig.row <- coef.match.name.df$panel == 'Gene-{MSigDB Collection}-Gene-Disease DWPC Feature'
coef.match.name.df[msig.row, 'easy_name'] <- paste('{', coef.match.name.df[msig.row, 'easy_name'], '}', sep='')


coef.match.name.df <- coef.match.name.df[match(gsub('XB|', '', XB.names, fixed=TRUE), coef.match.name.df$feature), ]


#prediction.df$gene_code <- hgnc.df[match(prediction.df$gene_symbol, hgnc.df$symbol), 'hgnc_id']

dir.create(file.path(webdata.dir, 'prediction-terms'), showWarnings=FALSE)
for (disease.code in gsub(':', '_', disease.summary.df$disease_code, fixed=TRUE)) {
  dir.create(file.path(webdata.dir, 'prediction-terms', disease.code), showWarnings=FALSE)
}

template.df <- data.frame('feature'=coef.match.name.df$easy_name)

for (i in 1:nrow(prediction.df)) {
  prediction.row <- prediction.df[i, ]
  disease.code <- gsub(':', '_', prediction.row$disease_code, fixed=TRUE)
  filename <- sprintf('%s.txt', gsub(':', '_', prediction.row$gene_code, fixed=TRUE))
  path <- file.path(webdata.dir, 'prediction-terms', disease.code, filename)
  prediction.terms <- cbind(template.df, 'value'=as.numeric(prediction.row[, XB.names]))
  prediction.terms <- format(prediction.terms, digits=2)
  write.table(prediction.terms, path, sep='\t', row.names=FALSE, quote=FALSE)
}


feature.description.df <- coef.match.name.df[, c('feature', 'name', 'metric', 'easy_name')]
path <- vtm.test.path <- file.path(modeling.dir, 'feature-descriptions.txt')
write.table(feature.description.df, path, sep='\t', row.names=FALSE, quote=FALSE)


