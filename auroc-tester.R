options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'

# Load code
code.dir <- file.path(project.dir, 'rcode')
sources <- c('settings.R', 'hnlp-learning.R', 'machine-learning.R', 'plotting.R')
for (source.filename in sources) {
  source(file.path(code.dir, source.filename))
}

# Load network info
network.id <- '140615-training'
network.dir <- file.path(project.dir, 'networks', network.id)

network.id <- '140615-0-training'
network.perm.dir <- file.path(project.dir, 'networks', 'permuted', network.id)

path <- file.path(network.dir, 'model', 'aucs.txt')
auc.df <- read.delim(path)
auc.df <- subset(auc.df, type == 'feature' & breadth == 'disease_specific')

path <- file.path(network.perm.dir, 'model', 'aucs.txt')
auc.perm.df <- read.delim(path)
auc.perm.df <- subset(auc.perm.df, type == 'feature' & breadth == 'disease_specific')

auc.df$auroc_perm <- auc.perm.df$auroc
auc.df$auprc_perm <- auc.perm.df$auprc

PairedT <- function(auc.df) {
  tt <- t.test(auc.df$auprc, auc.df$auprc_perm, alternative='greater', paired=FALSE)
  pval <- tt$p.value
  c('p_value'=pval, 'passed'=pval <= 0.05)
}

plyr::ddply(auc.df, 'feature', PairedT)




library(pAUC)


### Delong method
# Read Features
feature.df <- ReadFeatures(file.path(network.dir, 'features'))
feature.names <- colnames(feature.df)[-(1:8)]

feature.perm.df <- ReadFeatures(file.path(network.perm.dir, 'features'))













