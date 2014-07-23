options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'

# Load code
code.dir <- file.path(project.dir, 'rcode')
sources <- c('settings.R', 'hnlp-learning.R', 'machine-learning.R', 'plotting.R')
for (source.filename in sources) {
  source(file.path(code.dir, source.filename))
}

# Load network info
for (permutation.number in 0:4) {
network.date <- '140615'
network.id <- sprintf('%s-%s', network.date, permutation.number)
cat(paste0(network.id, '\n'))
network.dir <- file.path(project.dir, 'networks', 'permuted', network.id)
dirs <- InitializeNetworkDir(network.dir)


# Read Features
feature.df <- ReadFeatures(dirs$features)
feature.names <- colnames(feature.df)[-(1:9)]

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
fit.select <- TrainModel(X=Xe.train, y=y.train, alpha=0)

SaveFit(fit.ridge, dirs, suffix='-ridge')
SaveFit(fit.lasso, dirs, suffix='-lasso')
SaveFit(fit.select, dirs, suffix='-select')

# Create coefficient data.frame
ridge.coef.df <- GLMNetCoef(fit.ridge$cv.model, X.train, y.train, prepend='ridge_',
  name=c('intercept', feature.converter[feature.names]))
lasso.coef.df <- GLMNetCoef(fit.lasso$cv.model, X.train, y.train, prepend='lasso_')
select.coef.df <- GLMNetCoef(fit.select$cv.model, Xe.train, y.train, prepend='select_')
coef.df <- merge(ridge.coef.df, lasso.coef.df, select.coef.df)
coef.path <- file.path(dirs$model, 'coefficients.txt')
write.table(coef.df, coef.path, sep='\t', row.names=FALSE, quote=FALSE)

# Make predictions using the global model
feature.df[, 'ridge'] <- MakePredictions(cv.model=fit.ridge$cv.model, X=X)
feature.df[, 'lasso'] <- MakePredictions(cv.model=fit.lasso$cv.model, X=X)
feature.df[, 'select'] <- MakePredictions(cv.model=fit.select$cv.model, X=Xe)

# Calculate performance
vtm.ridge <- fit.ridge$vtm
vtm.lasso <- fit.lasso$vtm
vtm.select <- fit.select$vtm

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

path <- file.path(dirs$plots, 'AUROC.pdf')
OpenPDF(path, width=width.full, height=width.full)
PlotAUROCs(auroc.df)
ClosePDF(path)

}
