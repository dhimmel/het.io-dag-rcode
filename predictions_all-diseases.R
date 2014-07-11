options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'

# Load code
code.dir <- file.path(project.dir, 'rcode')
sources <- c('settings.R', 'hnlp-learning.R', 'machine-learning.R', 'plotting.R')
for (source.filename in sources) {
  source(file.path(code.dir, source.filename))
}

# Load network info
network.id <- '140615-all-assoc'
network.dir <- file.path(project.dir, 'networks', network.id)
dirs <- InitializeNetworkDir(network.dir)
dirs$webdata <- file.path(network.dir, 'webdata')
dir.create(dirs$webdata, showWarnings=FALSE)

# Read Features
feature.df <- ReadFeatures(dirs$features)
feature.names <- colnames(feature.df)[-(1:8)]

# Read feature descriptions and easy_names
feature.desc.path <- file.path(project.dir, 'data-integration', 'feature-descriptions.txt')
desc.df <- read.delim(feature.desc.path)
feature.converter <- desc.df$easy_name
names(feature.converter) <- desc.df$feature

# Read pathophysiology
pathophys.path <- file.path(project.dir, 'data-integration', 'pathophysiology.txt')
pathophys.df <- read.delim(pathophys.path)
feature.df$disease_pathophys <- pathophys.df[match(feature.df$disease_code, pathophys.df$disease_code), 'pathophysiology']

################################################################################
## Fit on all oberservations

# Fit on whole data (global)
X <- as.matrix(feature.df[, feature.names])
y <- feature.df$status_int

X.train <- as.matrix(subset(feature.df, status_int != -1, select=feature.names))
y.train <- subset(feature.df, status_int != -1)[, 'status_int']

fit.ridge <- TrainModel(X=X.train, y=y.train, alpha=0)
fit.lasso <- TrainModel(X=X.train, y=y.train, alpha=1)

SaveFit(fit.ridge, dirs, suffix='-ridge')
SaveFit(fit.lasso, dirs, suffix='-lasso')

# Create coefficient data.frame
ridge.coef.df <- GLMNetCoef(fit.ridge$cv.model, X.train, y.train, prepend='ridge_',
  name=c('intercept', feature.converter[feature.names]))
lasso.coef.df <- GLMNetCoef(fit.lasso$cv.model, X.train, y.train, prepend='lasso_')
coef.df <- merge(ridge.coef.df, lasso.coef.df)
coef.path <- file.path(dirs$model, 'coefficients.txt')
write.table(coef.df, coef.path, sep='\t', row.names=FALSE, quote=FALSE)

# Make predictions using the global model
feature.df[, 'ridge'] <- MakePredictions(cv.model=fit.ridge$cv.model, X=X)
feature.df[, 'lasso'] <- MakePredictions(cv.model=fit.lasso$cv.model, X=X)

# Calculate performance
vtm.ridge <- fit.ridge$vtm
vtm.lasso <- fit.lasso$vtm

# Save ridge predictions
prediction.cast <- reshape2::dcast(feature.df, gene_symbol ~ disease_name, value.var='ridge')
predictions.cast.path <- file.path(dirs$model, 'prediction-table.txt')
write.table(prediction.cast, predictions.cast.path, sep='\t', row.names=FALSE, quote=FALSE)

## Save XB and make predictions
XB.mat <- X %*% diag(ridge.coef.df[-1, 'ridge_coef'])
XB.names <- paste('XB|', feature.names, sep='')
colnames(XB.mat) <- XB.names
prediction.df <- cbind(feature.df, as.data.frame(XB.mat, check.names=FALSE))
predictions.file <- gzfile(file.path(dirs$model, 'predictions.txt.gz'), 'w')
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

path <- file.path(dirs$plots, 'coefficients.pdf')
OpenPDF(path, width=width.half, height=4.5)
print(gg.zcoef)
ClosePDF(path)


################################################################################
## Performance by status type

pos.statuses <- c('HC_primary', 'HC_secondary', 'LC_primary', 'LC_secondary')
auroc.list <- list()

roc.df <- do.call(rbind, lapply(pos.statuses, function(pos.status) {
  subset.df <- subset(feature.df, status %in% c('negative', pos.status))
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

path <- file.path(dirs$plots, 'ROC-by-positive-set.pdf')
OpenPDF(path, width=width.half, height=width.half - 0.1)
print(gg.roc)
ClosePDF(path)

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
cor.melt$label <- sprintf('%02d', round(cor.melt$correlation * 100))

xlimits <- features.sorted[-length(features.sorted)]
ylimits <- rev(features.sorted[-1])

diverging.cols <- c('#0000FF', '#f7eff7', '#FF0000') # blues
diverging.cols <- Solar('green', 'cyan', 'blue', 'violet', 'magenta', 'red')
cor.plot <- ggplot(cor.melt, aes(x, y, fill=correlation))
cor.plot <- SetGGTheme(cor.plot) +
  geom_tile(aes(alpha=abs(correlation)), color='white') + 
  geom_text(aes(label=label), size=3.4, color=Solar('base02')) + 
  scale_fill_gradientn(colours=diverging.cols) + 
  scale_alpha(range=c(0.3, 0.9)) + guides(alpha=FALSE, fill=FALSE) +
  scale_x_discrete(limits=xlimits, labels=feature.converter[xlimits]) + 
  scale_y_discrete(limits=ylimits, labels=feature.converter[ylimits]) +
  theme(axis.text.x=element_text(angle=65, hjust=1)) +
  xlab(NULL) + ylab(NULL) + coord_fixed() +
  theme(plot.background=element_blank(), panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(), panel.border=element_blank(),
    axis.line=element_line(color='black'))

path <- file.path(dirs$plots, 'feature-correlations.pdf')
OpenPDF(path, width=width.full, height=width.full)
print(cor.plot)
ClosePDF(path)

################################################################################
## Disease Specific Predictions on global data set

feat.perf.df <- subset(feature.df, status_int != -1)
auroc.df <- ComputeAUROCDF(feat.perf.df)
path <- file.path(dirs$model, 'aurocs.txt')
write.table(auroc.df, path, sep='\t', row.names=FALSE, quote=FALSE)
auroc.df <- read.delim(path)

path <- file.path(project.dir, 'networks', 'permuted', 'plots', 'AUROCs-permuted.txt')
auroc.perm.df <- read.delim(path)

path <- file.path(dirs$plots, 'AUROC.pdf')
OpenPDF(path, width=width.full, height=width.full)
PlotAUROCs(auroc.df, perm.df=auroc.perm.df)
ClosePDF(path)


################################################################################
# Webdata

#################
# Disease Summary Table
disease.summary.df <- subset(auroc.df, breadth == 'disease_specific' & name == 'Ridge',
  select=c('disease_name', 'disease_code', 'disease_pathophys', 'positives', 'auroc'))

colnames(disease.summary.df)[colnames(disease.summary.df) == 'positives'] <- 'associations'
colnames(disease.summary.df)[colnames(disease.summary.df) == 'disease_pathophys'] <- 'pathophysiology'
disease.summary.path <- file.path(dirs$webdata, 'disease-summary-table.txt')
disease.summary.df <- format(disease.summary.df, digits=2)
write.table(disease.summary.df, disease.summary.path, sep='\t', row.names=FALSE, quote=FALSE)

#################
# Gene Summary Table
gene.summary.df <- plyr::ddply(feature.df, c('gene_symbol', 'gene_code'), plyr::summarize,
  'associations'=sum(status_int == 1),
  'mean_prediction'=mean(ridge) * 100
)
gene.summary.path <- file.path(dirs$webdata, 'gene-summary-table.txt')
gene.summary.df <- format(gene.summary.df, digits=2)
write.table(gene.summary.df, gene.summary.path, sep='\t', row.names=FALSE, quote=FALSE)

#################
# Feature Summary Table
feature.summary.df <- subset(auroc.df, breadth == 'global' & type == 'feature',
  select=c('feature', 'name', 'metric', 'auroc'))
colnames(feature.summary.df)[colnames(feature.summary.df) == 'name'] <- 'metapath'
feature.summary.df[feature.summary.df$metric == 'PC_s', 'metric'] <- 'Path Count'
feature.summary.df[feature.summary.df$metric == 'PC_t', 'metric'] <- 'Path Count'
feature.summary.df[feature.summary.df$metric == 'DWPC_0.4', 'metric'] <- 'DWPC (w=0.4)'
feature.summary.df$standardized_coefficient <- coef.df[match(feature.summary.df$feature, coef.df$feature), 'ridge_zcoef']
feature.summary.df$feature <- gsub('|', '_', feature.summary.df$feature, fixed=TRUE)
feature.summary.path <- file.path(dirs$webdata, 'feature-summary-table.txt')
feature.summary.df <- format(feature.summary.df, digits=2)
write.table(feature.summary.df, feature.summary.path, sep='\t', row.names=FALSE, quote=FALSE)


#################
# Disease Tables

status.converter <- c('HC_primary'='+ HC-P', 'LC_primary'='\u00B1 LC-P', 'HC_secondary'='\u00B1 HC-S', 'LC_secondary'='\u00B1 LC-S', 'negative'='-')

dir.create(file.path(dirs$webdata, 'disease-tables'), showWarnings=FALSE)
MakeDiseaseTable <- function(feature.df) {
  gene.summary.index <- match(feature.df$gene_symbol, gene.summary.df$gene_symbol)
  disease.table <- data.frame(
    'gene_symbol'=feature.df$gene_symbol, 
    'gene_code'=feature.df$gene_code, 
    'status'=status.converter[feature.df$status], 
    'other_associations'=feature.df[, 'PC_s|G-a-D'],
    'mean_prediction'=gene.summary.df[gene.summary.index, 'mean_prediction'],
    'prediction'=feature.df$ridge * 100)
  return(disease.table)
}
disease.tables <- plyr::dlply(feature.df, 'disease_code', MakeDiseaseTable)
for (disease_code in names(disease.tables)) {
  disease.table <- disease.tables[[disease_code]]
  disease.table <- format(disease.table, digits=2)
  filename <- sprintf('%s.txt', gsub(':', '_', disease_code, fixed=TRUE))
  table.path <- file.path(dirs$webdata, 'disease-tables', filename)
  write.table(disease.table, table.path, sep='\t', row.names=FALSE, quote=FALSE)
}

#################
# Gene Tables
dir.create(file.path(dirs$webdata, 'gene-tables'), showWarnings=FALSE)
MakeGeneTable <- function(feature.df) {
  gene.table <- data.frame(
    'disease_name'=feature.df$disease_name, 
    'disease_code'=feature.df$disease_code, 
    'pathophysiology'=feature.df$disease_pathophys,
    'status'=status.converter[feature.df$status], 
    'other_associations'=feature.df[, 'PC_t|G-a-D'],
    'prediction'=feature.df$ridge * 100)
  return(gene.table)
}
gene.tables <- plyr::dlply(feature.df, 'gene_code', MakeGeneTable)
for (gene.code in names(gene.tables)) {
  gene.table <- gene.tables[[gene.code]]
  gene.table <- format(gene.table, digits=2)
  filename <- sprintf('%s.txt', gsub(':', '_', gene.code, fixed=TRUE))
  table.path <- file.path(dirs$webdata, 'gene-tables', filename)
  write.table(gene.table, table.path, sep='\t', row.names=FALSE, quote=FALSE)
}


#################
# Feature Tables
dir.create(file.path(dirs$webdata, 'feature-tables'), showWarnings=FALSE)
fXd.auroc.df <- subset(auroc.df, breadth == 'disease_specific' & type == 'feature')
fXd.auroc.df$feature <- gsub('|', '_', fXd.auroc.df$feature, fixed=TRUE)

MakeFeatureTable <- function(feature.auroc.df) {
  feature.table <- data.frame(
    'disease_name'=feature.auroc.df$disease_name,
    'disease_code'=feature.auroc.df$disease_code,
    'pathophysiology'=feature.auroc.df$disease_pathophys,
    'associations'=feature.auroc.df$positives,
    'model_auroc'=disease.summary.df[match(feature.auroc.df$disease_code, disease.summary.df$disease_code), 'auroc'],
    'auroc'=feature.auroc.df$auroc)
  return(feature.table)
}

feature.tables <- plyr::dlply(fXd.auroc.df, 'feature', MakeFeatureTable)
for (feature in names(feature.tables)) {
  feature.table <- feature.tables[[feature]]
  feature.table <- format(feature.table, digits=2)
  filename <- sprintf('%s.txt', feature)
  table.path <- file.path(dirs$webdata, 'feature-tables', filename)
  write.table(feature.table, table.path, sep='\t', row.names=FALSE, quote=FALSE)
}


#################
# Coef Tables

coef.match.name.df <- subset(auroc.df, breadth == 'global' & type == 'feature',
  select=c('feature', 'name', 'metric'))
coef.match.name.df <- coef.match.name.df[match(gsub('XB|', '', XB.names, fixed=TRUE), coef.match.name.df$feature), ]

dir.create(file.path(dirs$webdata, 'prediction-terms'), showWarnings=FALSE)
for (disease.code in gsub(':', '_', disease.summary.df$disease_code, fixed=TRUE)) {
  dir.create(file.path(dirs$webdata, 'prediction-terms', disease.code), showWarnings=FALSE)
}

template.df <- data.frame('feature'=coef.match.name.df$name)

for (i in 1:nrow(prediction.df)) {
  prediction.row <- prediction.df[i, ]
  disease.code <- gsub(':', '_', prediction.row$disease_code, fixed=TRUE)
  filename <- sprintf('%s.txt', gsub(':', '_', prediction.row$gene_code, fixed=TRUE))
  path <- file.path(dirs$webdata, 'prediction-terms', disease.code, filename)
  prediction.terms <- cbind(template.df, 'value'=as.numeric(prediction.row[, XB.names]))
  prediction.terms <- format(prediction.terms, digits=2)
  write.table(prediction.terms, path, sep='\t', row.names=FALSE, quote=FALSE)
}


