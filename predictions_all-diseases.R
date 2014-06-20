library(glmnet)
library(ggplot2)
library(doMC)
library(grid)
library(gtools)
library(plyr)
library(reshape2)

options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
network.id <- '140615-all-assoc'

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

# Read feature descriptions and easy_names
feature.desc.path <- file.path(project.dir, 'data-integration', 'feature-descriptions.txt')
desc.df <- read.delim(feature.desc.path)
feature.converter <- desc.df$easy_name
names(feature.converter) <- desc.df$feature

# Read pathophysiology
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

# Fit regression (with cross-validation to choose 'one-standard-error' lambda)
doMC::registerDoMC(cores=7)
set.seed(0); cv.ridge <- glmnet::cv.glmnet(X.model, y.model, family='binomial', 
  alpha=0, standardize=TRUE, parallel=TRUE) # fit ridge
set.seed(0); cv.lasso <- glmnet::cv.glmnet(X.model, y.model, family='binomial', 
  alpha=1, standardize=TRUE, parallel=TRUE) # fit lasso
# Save cv.glmnet objects as an R data structures
saveRDS(cv.ridge, file.path(modeling.dir, 'model-ridge.rds'))
saveRDS(cv.lasso, file.path(modeling.dir, 'model-lasso.rds'))

# Read ridge and lasso models from file
cv.ridge <- readRDS(file.path(modeling.dir, 'model-ridge.rds'))
cv.lasso <- readRDS(file.path(modeling.dir, 'model-lasso.rds'))

# Create coefficient data.frame
ridge.coef.df <- GLMNetCoef(cv.ridge, X.model, y.model, prepend='ridge_',
  name=c('intercept', feature.converter[feature.names]))
lasso.coef.df <- GLMNetCoef(cv.lasso, X.model, y.model, prepend='lasso_')
coef.df <- merge(ridge.coef.df, lasso.coef.df)
coef.path <- file.path(modeling.dir, 'coefficients.txt')
write.table(coef.df, coef.path, sep='\t', row.names=FALSE, quote=FALSE)


# Make predictions using the global model
features.df[, 'ridge'] <- as.numeric(predict(cv.ridge, s=cv.ridge$lambda.1se, newx=X, type='response'))
features.df[, 'lasso'] <- as.numeric(predict(cv.lasso, s=cv.lasso$lambda.1se, newx=X, type='response'))

# Calculate performance
feat.perf.df <- subset(features.df, status_int != -1)
vtm.ridge <- VariableThresholdMetrics(feat.perf.df$ridge, feat.perf.df$status_int)
vtm.lasso <- VariableThresholdMetrics(feat.perf.df$lasso, feat.perf.df$status_int)

# Save ridge predictions
prediction.cast <- reshape2::dcast(features.df, gene_symbol ~ disease_name, value.var='ridge')
predictions.cast.path <- file.path(modeling.dir, 'prediction-table.txt')
write.table(prediction.cast, predictions.cast.path, sep='\t', row.names=FALSE, quote=FALSE)

## Save XB and make predictions
XB.mat <- X %*% diag(ridge.coef.df[-1, 'ridge_coef'])
XB.names <- paste('XB|', feature.names, sep='')
colnames(XB.mat) <- XB.names
prediction.df <- cbind(features.df, as.data.frame(XB.mat, check.names=FALSE))
predictions.file <- gzfile(file.path(modeling.dir, 'predictions.txt.gz'), 'w')
write.table(prediction.df, predictions.file, sep='\t', row.names=FALSE, quote=FALSE)
close(predictions.file)


################################################################################
## Coefficient Plot

coef.df$neg_ridge_zcoef <- -coef.df$ridge_zcoef
coef.melt <- reshape2::melt(subset(coef.df, feature != 'intercept'),
  measure.vars=c('neg_ridge_zcoef', 'lasso_zcoef'))

coef.melt.ridge <- subset(coef.melt, variable=='neg_ridge_zcoef')
coef.sorted <- coef.melt.ridge[order(coef.melt.ridge$value, decreasing=TRUE), 'name']

gg.zcoef <- ggplot(coef.melt, aes(x=name, ymin=0, ymax=value, color=variable))
gg.zcoef <- SetGGTheme(gg.zcoef) +
  geom_linerange(stat='identity', size=4.2) + 
  geom_hline(yintercept=0, color=Solar('base02')) +
  coord_flip() +
  ylab('Standardized Coefficient') +
  scale_x_discrete(limits=coef.sorted) + 
  scale_y_continuous(labels=abs) +
  scale_color_manual(name='Method (AUROC)', values=c(Solar('red'), Solar('blue')), 
    labels=c(sprintf('ridge (%.3f)', vtm.ridge$auroc), 
             sprintf('lasso (%.3f)', vtm.lasso$auroc))) +
  theme(legend.key=element_rect(linetype='blank')) +
  theme(axis.text.y=element_text(angle=30, hjust=1)) + xlab(NULL) +
  theme(legend.justification=c(1,0), legend.position=c(1,0), 
    legend.background=element_rect(color='grey60', size=0.2))

ggsave(file.path(graphics.dir, 'coefficients.pdf'), gg.zcoef, width=width.half, height=4.5)


################################################################################
## Performance by status type

pos.statuses <- c('HC_primary', 'HC_secondary', 'LC_primary', 'LC_secondary')
auroc.list <- list()

roc.df <- do.call(rbind, lapply(pos.statuses, function(pos.status) {
  subset.df <- subset(features.df, status %in% c('negative', pos.status))
  status <- subset.df$status == pos.status
  vtm <- VariableThresholdMetrics(subset.df$ridge, status)
  roc.df <- vtm$roc.df
  roc.df[, 'Positives'] <- pos.status
  auroc.list[[pos.status]] <<- vtm$auroc
  return(roc.df)
}))


cols <- as.character(solarized[c('red', 'blue', 'red', 'blue')])
linetypes <- c('solid', 'solid', 'dotted', 'dotted')
gglabels <- c(
  sprintf('HC Primary (%.2f)', auroc.list$HC_primary),
  sprintf('HC Secondary (%.2f)', auroc.list$HC_secondary),
  sprintf('LC Primary (%.2f)', auroc.list$LC_primary),
  sprintf('LC Secondary (%.2f)', auroc.list$LC_secondary))

gg.roc <- ggplot(roc.df, aes(fpr, recall, color=Positives, linetype=Positives))
gg.roc <- ggROC(gg.roc) + geom_path(size=0.9) +
  theme(legend.title=element_text(size=9.5), legend.text=element_text(size=9)) +
  scale_linetype_manual(name='Positive Set (AUROC)', 
    values=linetypes, labels=gglabels, breaks=pos.statuses) +
  scale_color_manual(name='Positive Set (AUROC)', 
    values=cols, labels=gglabels, breaks=pos.statuses)
ggsave(file.path(graphics.dir, 'ROC-by-positive-set.pdf'), gg.roc, width=width.half, height=width.half - 0.1)


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
auroc.df1 <- rbind(
plyr::ddply(feat.perf.df, c('disease_code'), plyr::summarize,
  'breadth' = 'disease_specific',
  'type'='model', 'name'='Ridge',
  'auroc' = VariableThresholdMetrics(ridge, status_int)$auroc,
  'positives'=sum(status_int)),
plyr::ddply(feat.perf.df, c('disease_code'), plyr::summarize,
  'breadth' = 'disease_specific',
  'type'='model', 'name'='Lasso',
  'auroc' = VariableThresholdMetrics(lasso, status_int)$auroc,
  'positives'=sum(status_int))
)
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


# Global AUC for Ridge Model
auroc.df4 <- rbind(
data.frame('breadth'='global', 'type'='model', 'name'='Ridge',
  'positives'=sum(feat.perf.df$status_int), 'auroc'=vtm.ridge$auroc),
data.frame('breadth'='global', 'type'='model', 'name'='Lasso',
  'positives'=sum(feat.perf.df$status_int), 'auroc'=vtm.lasso$auroc)
)
auroc.df4 <- OrderAurocDf(auroc.df4)

# Add positives to auroc.df3
auroc.df3$positives <- auroc.df4[match(auroc.df3$disease_code, auroc.df4$disease_code), 'positives']
auroc.df3 <- OrderAurocDf(auroc.df3)

# Combine
auroc.df <- rbind(auroc.df1, auroc.df2, auroc.df3, auroc.df4)
auroc.df$disease_pathophys <- pathophys.df[match(auroc.df$disease_code, pathophys.df$disease_code), 'pathophysiology']
auroc.df$disease_name <- pathophys.df[match(auroc.df$disease_code, pathophys.df$disease_code), 'disease_name']

metric.metapath.mat <- do.call(rbind, strsplit(auroc.df$feature, '|', fixed=TRUE))
auroc.df$metric <- metric.metapath.mat[, 1]
auroc.df$metapath <- gsub('-', '', metric.metapath.mat[, 2])
auroc.df[is.na(auroc.df$name), 'name'] <- feature.converter[auroc.df[is.na(auroc.df$name), 'feature']]


################################################################################
## Plot AUROCs

NAtoFALSE <- function(vec) {
  vec[is.na(vec)] <- FALSE
  return(vec)
}
auroc.df$panel <- 'Model'
auroc.df[NAtoFALSE(substr(auroc.df$metric, 1, 4)) == 'PC_s', 'panel'] <- 'PCs'
auroc.df[NAtoFALSE(substr(auroc.df$metric, 1, 4)) == 'PC_t', 'panel'] <- 'PCt'
auroc.df[NAtoFALSE(substr(auroc.df$metric, 1, 4)) == 'DWPC', 'panel'] <- 'DWPC'
is.msig.auc <- NAtoFALSE(substr(auroc.df$name, 1, 1) == '{')
auroc.df[is.msig.auc, 'panel'] <- 'Gene-{MSigDB Collection}-Gene-Disease DWPC Feature'

MeanConfInt <- function(x) {t.test(x)$conf.int[1:2]}
pathophys.colors <- c('#005200', '#B20000', '#8F008F', '#0000B2', '#E68A00', 'black')

#pathophys.colors <- as.character(solarized[c('red', 'yellow', 'magenta', 'violet', 'cyan', 'base01')])
ymin.auroc <- min(subset(auroc.df, auroc != 0)$auroc)


PerfPlot <- function(gg) {
  gg <- SetGGTheme(gg) +
  facet_grid(. ~ panel, scales='free_x', space='free_x') +
  geom_hline(aes(yintercept=0.5), color='grey', linetype='solid') +
  stat_summary(fun.y='MeanConfInt', geom='line', color='grey', size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  xlab(NULL) + ylab('AUROC') +
  scale_y_continuous(limits=c(ymin.auroc, 1), breaks=seq(0, 1, .2), expand=c(0.03, 0)) +
  scale_colour_manual(values=pathophys.colors, name='Pathophysiology')
  return(gg)
}

hyphen.size <- 15.5
jitter.width <- 0.3
point.size <- 1.75

## MSigDB Plot
msig.df <- subset(auroc.df, panel == 'Gene-{MSigDB Collection}-Gene-Disease DWPC Feature' | panel == 'Model')
msig.df$name <- gsub('[{}]', '', msig.df$name)
msig.disease.df <- subset(msig.df, breadth == 'disease_specific')
msig.global.df <- subset(msig.df, breadth == 'global')
# Order by global AUROC
msig.levels <- as.character(msig.global.df$name[order(msig.global.df$auroc)])
msig.disease.df$name <- factor(msig.disease.df$name, levels=msig.levels)
msig.global.df$name <- factor(msig.global.df$name, levels=msig.levels)

set.seed(0); msig.plot <- ggplot(msig.disease.df, aes(name, auroc)) 
msig.plot <- PerfPlot(msig.plot) +
  geom_point(data=msig.global.df, size=hyphen.size, shape='-', color=Solar('base02')) + 
  geom_point(aes(color=disease_pathophys), position=position_jitter(width=jitter.width), 
    alpha=0.7, size=point.size, show_guide=FALSE)

## NonMSigDB Plot
nonmsig.df <- subset(auroc.df, panel != 'Gene-{MSigDB Collection}-Gene-Disease DWPC Feature')
nonmsig.df$panel <- factor(nonmsig.df$panel, levels=c('DWPC', 'PCt', 'PCs', 'Model'))
nonmsig.disease.df <- subset(nonmsig.df, breadth == 'disease_specific' & panel != 'PCt')
nonmsig.global.df <- subset(nonmsig.df, breadth == 'global')
# Order by global AUROC
nonmsig.levels <- unique(as.character(nonmsig.global.df$name[order(nonmsig.global.df$auroc)]))
nonmsig.disease.df$name <- factor(nonmsig.disease.df$name, levels=nonmsig.levels)
nonmsig.global.df$name <- factor(nonmsig.global.df$name, levels=nonmsig.levels)

set.seed(0); nonmsig.plot <- ggplot(nonmsig.disease.df, aes(name, auroc))
nonmsig.plot <- PerfPlot(nonmsig.plot) +
  geom_point(data=nonmsig.global.df, size=hyphen.size, shape='-', color=Solar('base02')) + 
  geom_point(aes(color=disease_pathophys), position=position_jitter(width=jitter.width),
    alpha=0.7, size=point.size) + 
  theme(legend.key=element_rect(linetype='blank')) +
  theme(legend.margin=grid::unit(1, 'points'))

pdf(file.path(graphics.dir, 'AUROC.pdf'), width=width.full, height=width.full)
gridExtra::grid.arrange(msig.plot, nonmsig.plot, ncol=1, widths=c(1, 1))
dev.off()


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

status.converter <- c('HC_primary'='+ HC-P', 'LC_primary'='\u00B1 LC-P', 'HC_secondary'='\u00B1 HC-S', 'LC_secondary'='\u00B1 LC-S', 'negative'='-')

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


