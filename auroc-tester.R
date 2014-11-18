library(pROC)
library(doMC)

options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'

# Load code
code.dir <- file.path(project.dir, 'rcode')
sources <- c('settings.R', 'hnlp-learning.R', 'machine-learning.R', 'plotting.R')
for (source.filename in sources) {
  source(file.path(code.dir, source.filename))
}


##################################################################################################
## Training

# Load network info
network.id <- '140615-training'
network.dir <- file.path(project.dir, 'networks', network.id)

network.id <- '140615-0-training'
network.perm.dir <- file.path(project.dir, 'networks', 'permuted', network.id)

# Read Features
feature.df <- ReadFeatures(file.path(network.dir, 'features'))
feature.df <- subset(feature.df, part=='train' & status_int != -1)
feature.names <- colnames(feature.df)[-(1:8)]

feature.perm.df <- ReadFeatures(file.path(network.perm.dir, 'features'))
feature.perm.df <- subset(feature.perm.df, part=='train' & (status == 'negative' | network_status))

# Compare permuted and unpermuted feature AUROCs

CompareROCs <- function(feature) {
  # pROC Documentation: doi:10.1186/1471-2105-12-77
  # http://cran.at.r-project.org/web/packages/pROC/pROC.pdf
  cat(paste0('Comparing ROCs for ', feature, '\n'))

  roc1 <- pROC::roc(response=feature.df$status_int, predictor=feature.df[, feature])
  roc0 <- pROC::roc(response=feature.perm.df$network_status, predictor=feature.perm.df[, feature])
  test <- pROC::roc.test(roc1=roc1, roc2=roc0, 
    method='delong', paired=FALSE, alternative='greater')

  test.df <- data.frame(
    'feature'=feature, 'name'=feature.converter[feature],
    'auroc'=roc1$auc, 'auroc_perm'=roc0$auc, 'p_value'=test$p.value,
    row.names=feature)
  return(test.df)
}


doMC::registerDoMC(cores=7)
roctest.df <- plyr::adply(feature.names, 1, CompareROCs, .parallel=TRUE)
roctest.df$X1 <- NULL # remove column
roctest.df$select <- roctest.df$p_value <= 0.05 | roctest.df$feature %in% c('PC_s|G-a-D', 'PC_t|G-a-D')
#roctest.df$name <- as.character(feature.converter[feature.names])

for (directory in c(network.dir, network.perm.dir)) {
  path <- file.path(directory, 'model', 'enhancing-features.txt')
  write.table(roctest.df, path, sep='\t', row.names=FALSE, quote=FALSE)
}

##################################################################################################

# Testing ROC comparison for predictions
path <- file.path(network.dir, 'model', 'testing-predictions-ridge.txt.gz')
feature.df <- read.delim(path, check.names=FALSE)

path <- file.path(network.perm.dir, 'model', 'testing-predictions-ridge.txt.gz')
feature.perm.df <- read.delim(path, check.names=FALSE)

feature.df$status_int <- feature.df$status
feature.perm.df$network_status <- feature.perm.df$status

testing.roc.df <- CompareROCs('prediction')

##################################################################################################
## All-assoc
network.id <- '140615-all-assoc'
network.dir <- file.path(project.dir, 'networks', network.id)

network.id <- '140615-0'
network.perm.dir <- file.path(project.dir, 'networks', 'permuted', network.id)

# Read Features
feature.df <- ReadFeatures(file.path(network.dir, 'features'))
feature.df <- subset(feature.df, status_int != -1)
feature.names <- colnames(feature.df)[-(1:8)]

feature.perm.df <- ReadFeatures(file.path(network.perm.dir, 'features'))

roctest.df <- plyr::adply(feature.names, 1, CompareROCs, .parallel=TRUE)
roctest.df$X1 <- NULL # remove column
roctest.df$select <- roctest.df$p_value <= 0.05 | roctest.df$feature %in% c('PC_s|G-a-D', 'PC_t|G-a-D')
#roctest.df$name <- as.character(feature.converter[feature.names])

for (directory in c(network.dir, network.perm.dir)) {
  path <- file.path(directory, 'model', 'enhancing-features.txt')
  write.table(roctest.df, path, sep='\t', row.names=FALSE, quote=FALSE)
}

