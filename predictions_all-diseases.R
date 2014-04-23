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
code.dir <- file.path(project.dir, 'rcode')
network.dir <- file.path(project.dir, 'networks', '140321-all-assoc')
feature.dir <- file.path(network.dir, 'features')
graphics.dir <- file.path(network.dir, 'graphics')
modeling.dir <- file.path(network.dir, 'modeling')
webdata.dir <- file.path(network.dir, 'webdata')
dir.create(graphics.dir, showWarnings=FALSE)
dir.create(modeling.dir, showWarnings=FALSE)
dir.create(webdata.dir, showWarnings=FALSE)

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
  colnames(feature.df)[colnames(feature.df) == 'source'] <- 'gene_symbol'
  colnames(feature.df)[colnames(feature.df) == 'target'] <- 'disease_code'
  colnames(feature.df)[colnames(feature.df) == 'target_name'] <- 'disease_name'
  feature.df.list[[doid_code]] <- feature.df
}

features.df <- do.call(rbind, feature.df.list)
feature.names <- colnames(feature.df)[-(1:5)]
category.path <- file.path(project.dir, 'data-integration', 'disease-categories.txt')
category.df <- read.delim(category.path)
features.df$disease_category <- category.df[match(features.df$disease_code, category.df$doid_code), 'category']


################################################################################
## Training and Testing Partitioning

train.df <- subset(features.df, part == 'train')
test.df <- subset(features.df, part == 'test')

X.train <- as.matrix(train.df[, feature.names])
X.test <- as.matrix(test.df[, feature.names])
y.train <- train.df$status
y.test <- test.df$status

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


# Testing Precision-Recall Curve
prc.plot <- ggplot(vtm.test.df[nrow(vtm.part.df):1, ], aes(recall, precision, color=threshold)) + 
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
roc.plot <- ggplot(vtm.part.df, aes(fpr, recall, linetype=part)) + 
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

################################################################################
## Refit on all oberservations

# Fit on whole data (global)
X <- as.matrix(features.df[, feature.names])
y <- features.df$status
set.seed(0)
cv.ridge <- glmnet::cv.glmnet(X, y, family='binomial', 
  alpha=glmnet.alpha, standardize=TRUE, parallel=TRUE)
# Save cv.glmnet object as an R data structure
saveRDS(cv.ridge, file.path(modeling.dir, 'global-model.rds'))
#cv.ridge <- readRDS(file.path(modeling.dir, 'global-model.rds'))
lambda.global <- cv.ridge$lambda.1se

# Make predictions using the global model
## Save coefficients
coefs.global <- coef(cv.ridge, s=lambda.global)[, 1]
zcoefs.global <- c(0, coefs.global[-1] * apply(X, 2, sd) / sd(y))
feature.name.mat <- do.call(rbind, strsplit(feature.names, '[|]'))
coef.df <- data.frame(
  'feature'=c('intercept', feature.names), 
  'metric'=c(NA, feature.name.mat[, 1]), 
  'metapath'=c(NA, feature.name.mat[, 2]),
  'coef'=coefs.global,
  'zcoef'=zcoefs.global, row.names=NULL)
coef.path <- file.path(modeling.dir, 'coefs-global.txt')
write.table(coef.df, coef.path, sep='\t', row.names=FALSE, quote=FALSE)

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
vtm.global <- VariableThresholdMetrics(features.df$prediction, features.df$status)

prediction.cast <- reshape2::dcast(features.df, gene_symbol ~ disease_name, value.var='prediction')
predictions.cast.path <- file.path(modeling.dir, 'prediction-table.txt')
write.table(prediction.cast, predictions.cast.path, sep='\t', row.names=FALSE, quote=FALSE)

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
auroc.df1 <- plyr::ddply(features.df, c('disease_code'), plyr::summarize,
  'breadth' = 'disease_specific',
  'type'='model', 'name'='ridge',
  'auroc' = VariableThresholdMetrics(prediction, status)$auroc,
  'positives'=sum(status))
auroc.df1 <- OrderAurocDf(auroc.df1)

ComputeFeatureAUC <- function(feat.df) {
  apply(feat.df[, feature.names], 2, function(x) VariableThresholdMetrics(x, feat.df$status)$auroc)
}

# Calculate global AUCs for each feature
auroc.vec2 <- ComputeFeatureAUC(features.df)
auroc.df2 <- data.frame(
  'feature'=names(auroc.vec2), 
  'auroc'=as.numeric(auroc.vec2),
  'type'='feature',
  'breadth'='global',
  'positives'=sum(features.df$status)
)
auroc.df2 <- OrderAurocDf(auroc.df2)

# Calculate disease-specific AUROCs for each feature
fXd.auroc.tab <- plyr::ddply(features.df, c('disease_code'), ComputeFeatureAUC)
auroc.df3 <- reshape2::melt(fXd.auroc.tab, id.vars=c('disease_code'), variable.name='feature', value.name='auroc')
auroc.df3 <- cbind(auroc.df3, 'type'='feature', 'breadth'='disease_specific')
auroc.df3$positives <- category.df[match(auroc.df3$disease_code, category.df$doid_code), 'count']
auroc.df3 <- OrderAurocDf(auroc.df3)

# Global AUC for Ridge Model
auroc.df4 <- data.frame('breadth'='global', 'type'='model', 'name'='ridge',
  'positives'=sum(features.df$status), 'auroc'=vtm.global$auroc)
auroc.df4 <- OrderAurocDf(auroc.df4)

# Combine
auroc.df <- rbind(auroc.df1, auroc.df2, auroc.df3, auroc.df4)
auroc.df$disease_name <- category.df[match(auroc.df$disease_code, category.df$doid_code), 'doid_name']
auroc.df$disease_category <- category.df[match(auroc.df$disease_code, category.df$doid_code), 'category']
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

## MSigDB Plot
msig.df <- subset(auroc.df, panel == 'Gene-{MSigDB Collection}-Gene-Disease DWPC Feature')# | panel == 'Model')
msig.disease.df <- subset(msig.df, breadth == 'disease_specific')
msig.global.df <- subset(msig.df, breadth == 'global')
# Order by global AUROC
msig.levels <- as.character(msig.global.df$name[order(msig.global.df$auroc)])
msig.disease.df$name <- factor(msig.disease.df$name, levels=msig.levels)
msig.global.df$name <- factor(msig.global.df$name, levels=msig.levels)

set.seed(0); msig.plot <- ggplot(msig.disease.df, aes(name, auroc)) + 
  facet_grid(. ~ panel, scales='free_x', space='free_x') +
  geom_hline(aes(yintercept=0.5), color='red', linetype='dashed') +
  stat_summary(fun.y='MeanConfInt', geom='line', color='grey', size=9.5) +
  geom_point(data=msig.global.df, size=11.5, shape='-', color='#4C4C4C') + 
  geom_point(aes(color=disease_category), position=position_jitter(width=0.28), alpha=0.7, size=2.5) + 
  theme(plot.margin = grid::unit(c(0, 0, 0, 0), 'points')) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  xlab(NULL) + ylab('AUROC') +# guides(color=FALSE) +
  theme(legend.position='top', legend.margin=grid::unit(-10, 'points')) +
  #theme(legend.justification=c(0,1), legend.position=c(0, 1),
  #  legend.background=element_rect(color='grey60', size=0.25),
  #  legend.margin=grid::unit(1, 'points'), legend.direction='horizontal') +
  scale_colour_manual(name='Disease Category', values=category.colors)
msig.plot.path <- file.path(graphics.dir, 'AUROC-msig-plot.pdf')
ggsave(msig.plot.path, msig.plot, width=6.83, height=4)

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
  geom_hline(aes(yintercept=0.5), color='red', linetype='dashed') +
  stat_summary(fun.y='MeanConfInt', geom='line', color='grey', size=12) +
  geom_point(data=nonmsig.global.df, size=15, shape='-', color='#4C4C4C') + 
  geom_point(aes(color=disease_category), position=position_jitter(width=0.3), alpha=0.7, size=2.5) + 
  theme(plot.margin = grid::unit(c(0, 0, 0, 0), 'points')) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  xlab(NULL) + ylab('AUROC') + 
  theme(legend.position='top', legend.margin=grid::unit(-10, 'points')) +
  #theme(legend.justification=c(0,1), legend.position=c(0, 1),
  #  legend.background=element_rect(color='grey60', size=0.25),
  #  legend.margin=grid::unit(1, 'points'), legend.direction='horizontal') +
  scale_colour_manual(name='Disease Category', values=category.colors)
nonmsig.plot.path <- file.path(graphics.dir, 'AUROC-nomsig-plot.pdf')
ggsave(nonmsig.plot.path, nonmsig.plot, width=6.83, height=4)


################################################################################
# Webdata
#https://stackoverflow.com/questions/5448545/how-to-retrieve-get-parameters-from-javascript
#http://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=1688
#https://www.datatables.net/forums/discussion/5611/how-to-grab-datatables-data-from-a-google-spreadsheet/p1

# https://code.google.com/p/jquery-csv/
# $.csv.toArray()
# http://datatables.net/ref#aaData

hgnc.path <- '/home/dhimmels/Documents/serg/data-sources/hgnc/140205/protein-coding.txt'
hgnc.df <- read.delim(hgnc.path)
features.df$gene_code <- hgnc.df[match(features.df$gene_symbol, hgnc.df$symbol), 'hgnc_id']

#################
# Disease Summary Table
# disease_name, disease_code, disease_category, associations, auroc
disease.summary.df <- subset(auroc.df, breadth == 'disease_specific' & type == 'model',
  select=c('disease_name', 'disease_code', 'disease_category', 'positives', 'auroc'))
#plyr::ddply(features.df, c('disease_code'), plyr::summarize, mean_prediction=mean(prediction))
#disease.summary.df$mean_prediction <- NA
colnames(disease.summary.df)[colnames(disease.summary.df) == 'positives'] <- 'associations'
disease.summary.path <- file.path(webdata.dir, 'disease-summary-table.txt')
disease.summary.df <- format(disease.summary.df, digits=2)
write.table(disease.summary.df, disease.summary.path, sep='\t', row.names=FALSE, quote=FALSE)

#################
# Gene Summary Table
#gene_symbol, gene_code, associations, mean_prediction
gene.summary.df <- plyr::ddply(features.df, c('gene_symbol'), plyr::summarize,
  'gene_code' = NA,
  'associations'=sum(status),
  'mean_prediction'=mean(prediction)
)
gene.summary.df$gene_code <- hgnc.df[match(gene.summary.df$gene_symbol, hgnc.df$symbol), 'hgnc_id']
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
feature.summary.df$standardized_coefficient <- coef.df[match(feature.summary.df$feature, coef.df$feature), 'zcoef']
feature.summary.df$feature <- gsub('|', '_', feature.summary.df$feature, fixed=TRUE)
feature.summary.path <- file.path(webdata.dir, 'feature-summary-table.txt')
feature.summary.df <- format(feature.summary.df, digits=2)
write.table(feature.summary.df, feature.summary.path, sep='\t', row.names=FALSE, quote=FALSE)


#################
# Disease Tables
dir.create(file.path(webdata.dir, 'disease-tables'), showWarnings=FALSE)
#gene_symbol, gene_code, positives, mean_prediction; prediction
MakeDiseaseTable <- function(feature.df) {
  gene.summary.index <- match(feature.df$gene_symbol, gene.summary.df$gene_symbol)
  disease.table <- data.frame(
    'gene_symbol'=feature.df$gene_symbol, 
    'gene_code'=gene.summary.df[gene.summary.index, 'gene_code'], 
    'status'=feature.df$status, 
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
    'disease_category'=feature.df$disease_category,
    'status'=feature.df$status, 
    'other_associations'=feature.df[, 'PC_t|G-a-D'],
    #'mean_prediction'=disease.summary.df[disease.summary.index, 'mean_prediction'],
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
    'disease_category'=feature.auroc.df$disease_category,
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


prediction.df$gene_code <- hgnc.df[match(prediction.df$gene_symbol, hgnc.df$symbol), 'hgnc_id']

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


