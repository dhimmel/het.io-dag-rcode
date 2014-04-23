library(ggplot2)
library(gridExtra)

directory <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/140313-thresholdsweep'
aucs.path <- file.path(directory, 'feature-aucs.txt')
auc.df <- read.delim(aucs.path)

threshold.DsD <- 0.2
threshold.GfG <- 0.2
threshold.GeT <- 1
threshold.DcT <- 32

# semantic similarity
df.DsD <- subset(auc.df, metapath %in% c('GaDsD', 'GaDsDsD'))
plot.DsD <- ggplot(df.DsD, aes(DsD, auc_logreg, shape=metapath)) +
  geom_vline(xintercept=threshold.DsD, color='darkgreen', linetype='dashed') + 
  geom_point() + geom_line() + theme_bw() + ylab('AUC') +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

# function
df.GfG <- subset(auc.df, metapath %in% c('GfGaD', 'GfGfGaD'))
plot.GfG <- ggplot(df.GfG, aes(GfG, auc_logreg, shape=metapath)) +
  geom_vline(xintercept=threshold.GfG, color='darkgreen', linetype='dashed') + 
  geom_point() + geom_line() + theme_bw() + ylab('AUC') +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

# expression
df.GeT <- subset(auc.df, metapath %in% c('GeTeGaD'))
plot.GeT <- ggplot(df.GeT, aes(GeT, auc_logreg, shape=metapath)) +
  geom_vline(xintercept=threshold.GeT, color='darkgreen', linetype='dashed') + 
  geom_point() + geom_line() + theme_bw() + ylab('AUC') +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

# cooccurrence 
df.DcT <- subset(auc.df, metapath %in% c('GaDcTcD'))
plot.DcT <- ggplot(df.DcT, aes(DcT, auc_logreg, shape=metapath)) +
  geom_vline(xintercept=threshold.DcT, color='darkgreen', linetype='dashed') + 
  geom_point() + geom_line() + theme_bw() + ylab('AUC') +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

# tissue localization
df.GeTcD <- subset(auc.df, metapath %in% c('GeTcD'))
plot.GeTcD <- ggplot(df.GeTcD, aes(GeT, TcD, fill=auc_global)) +
  geom_tile() + scale_fill_gradient(low='white', high='darkblue') +
  geom_vline(xintercept=threshold.GeT, color='darkgreen', linetype='dashed') + 
  geom_hline(yintercept=threshold.DcT, color='darkgreen', linetype='dashed') + 
  theme_bw() + guides(fill=guide_colorbar(title='GeTcD AUC')) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

pdf(file.path(directory, 'AUC-lineplots.pdf'), width=4, height=7)
gridExtra::grid.arrange(plot.DsD, plot.GeT, plot.DcT, ncol=1)
dev.off()

pdf(file.path(directory, 'AUC-tileplot.pdf'), width=5, height=4)
plot.GeTcD
dev.off()

#gridExtra::grid.arrange(plot.GeT, plot.GeTcD, ncol=1)

