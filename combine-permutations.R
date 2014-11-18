options(stringsAsFactors=FALSE)

dirs <- list()
dirs$project <- project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
dirs$permnets <- file.path(dirs$project, 'networks', 'permuted')
dirs$plots <- file.path(dirs$permnets, 'plots')
for (directory in c(dirs$plots)) {
  dir.create(directory, showWarnings=FALSE)
}

# Load code
code.dir <- file.path(dirs$project, 'rcode')
sources <- c('settings.R', 'machine-learning.R', 'plotting.R')
for (source.filename in sources) {
  print(source.filename)
  source(file.path(code.dir, source.filename))
}

# Load network info
network.date <- '140615'
perm.nums <- 0:4
network.ids <- NULL

auc.dfs <- list()
for (perm.num in perm.nums) {
  network.id <- sprintf('%s-%s', network.date, perm.num)
  network.ids <- c(network.ids, network.id)
  path <- file.path(dirs$permnets, network.id, 'model', 'aucs.txt')
  auc.dfs[[network.id]] <- read.delim(path)
}

mean.df <- auc.dfs[[1]]
mean.df$auroc <- apply(sapply(auc.dfs, function(auc.df) auc.df$auroc), 1, mean)
mean.df$auprc <- apply(sapply(auc.dfs, function(auc.df) auc.df$auprc), 1, mean)


path <- file.path(dirs$plots, 'AUCs-permuted.txt')
write.table(mean.df, path, sep='\t', row.names=FALSE, quote=FALSE)

path <- file.path(dirs$plots, 'AUROC-permuted-average.pdf')
OpenPDF(path, width=width.full, height=width.full)
PlotAUROCs(mean.df)
ClosePDF(path)

