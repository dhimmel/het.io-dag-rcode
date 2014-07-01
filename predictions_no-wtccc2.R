options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'

# Load code
code.dir <- file.path(project.dir, 'rcode')
sources <- c('settings.R', 'hnlp-learning.R', 'machine-learning.R', 'plotting.R')
for (source.filename in sources) {
  source(file.path(code.dir, source.filename))
}

# Load network info
network.id <- '140615-no-wtccc2'
network.dir <- file.path(project.dir, 'networks', network.id)
dirs <- InitializeNetworkDir(network.dir)

# Read Features
feature.df <- ReadFeatures(dirs$features)
feature.names <- colnames(feature.df)[-(1:8)]

# Fit on whole data (global)
X <- as.matrix(feature.df[, feature.names])
y <- feature.df$status_int

X.train <- as.matrix(subset(feature.df, status_int != -1, select=feature.names))
y.train <- subset(feature.df, status_int != -1)[, 'status_int']

fit.ridge <- TrainModel(X=X.train, y=y.train, alpha=0)
SaveFit(fit.ridge, dirs)

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
ms.novel.df <- subset(ms.df, status_int != 1 | gene_symbol %in% genes.wtccc2.novel)

X.test <- as.matrix(subset(ms.novel.df, status_int != -1 | wtc_novel, select=feature.names))
y.test <- subset(ms.novel.df, status_int != -1 | wtc_novel)[, 'wtc_novel']

test.ridge <- TestModel(cv.model=fit.ridge$cv.model, X=X.test, y=y.test)
SaveTest(test.ridge, dirs)


# Training and Testing ROC Curves
gg.roc <- ggplot(test.ridge$vtm$roc.df, aes(fpr, recall, color=part))
gg.roc <- ggROC(gg.roc) + geom_path(color=as.character(solarized['violet'])) +
  annotate('text', x=Inf, y=-Inf, 
    label=sprintf('AUROC == %.3f', test.ridge$vtm$auroc), 
    hjust=1.1, vjust=-1, parse=TRUE, size=3.5)

# Testing Precision-Recall Curve
prc.df <- PrunePRC(test.ridge$vtm$threshold.df)
gg.prc <- ggplot(prc.df[nrow(prc.df):1, ], aes(recall, precision, color=threshold))
ylim.prc <- range(prc.df$precision)
gg.prc <- ggPRC(gg.prc) +
  annotate('text', x=0.7, y=Inf, 
    label=sprintf('AUPRC == %s', ChrRound(test.ridge$vtm$auprc, 3)), 
    hjust=1, vjust=2, parse=TRUE, size=3.5)

# Save combined ROC and PRC
path <- file.path(dirs$plots, 'performance.pdf')
OpenPDF(path, width=width.full, height=2.5)
gridExtra::grid.arrange(gg.roc, gg.prc, nrow=1, widths=c(1, 1.625))
ClosePDF(path)

# MS Predictions for all genes
ms.df$ridge <- MakePredictions(cv.model=fit.ridge$cv.model,
  X=as.matrix(ms.df[, feature.names]))
path <- file.path(dirs$model, 'MS-predictions.txt.gz')
gz.file <- gzfile(path, 'w')
write.table(ms.df, gz.file, sep='\t', row.names=FALSE, quote=FALSE)
close(gz.file)

