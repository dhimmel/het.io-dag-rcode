library(ggplot2)
library(gridExtra)

directory <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/140518-thresholdsweep'
aucs.path <- file.path(directory, 'feature-aucs.txt')
auc.df <- read.delim(aucs.path)

#threshold.DsD <- 0.2
#threshold.GfG <- 0.2
threshold.GeT <- 1.4
threshold.DcT <- 29

# semantic similarity
#df.DsD <- subset(auc.df, metapath %in% c('GaDsD', 'GaDsDsD'))
#plot.DsD <- ggplot(df.DsD, aes(DsD, auc_logreg, shape=metapath)) +
#  geom_vline(xintercept=threshold.DsD, color='darkgreen', linetype='dashed') + 
#  geom_point() + geom_line() + theme_bw() + ylab('AUC') +
#  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

# function
#df.GfG <- subset(auc.df, metapath %in% c('GfGaD', 'GfGfGaD'))
#plot.GfG <- ggplot(df.GfG, aes(GfG, auc_logreg, shape=metapath)) +
#  geom_vline(xintercept=threshold.GfG, color='darkgreen', linetype='dashed') + 
#  geom_point() + geom_line() + theme_bw() + ylab('AUC') +
#  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

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

pdf(file.path(directory, 'AUC-lineplots-ridge.pdf'), width=4, height=4)
gridExtra::grid.arrange(plot.GeT, plot.DcT, ncol=1)
dev.off()


# tissue localization
df.GeTcD <- subset(auc.df, metapath %in% c('GeTcD'))
plot.GeTcD <- ggplot(df.GeTcD, aes(GeT, TcD, fill=auc_global)) +
  geom_tile() + scale_fill_gradient(low='white', high='darkblue') +
  geom_vline(xintercept=threshold.GeT, color='darkgreen', linetype='dashed') + 
  geom_hline(yintercept=threshold.DcT, color='darkgreen', linetype='dashed') + 
  theme_bw() + guides(fill=guide_colorbar(title='GeTcD AUC')) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

pdf(file.path(directory, 'AUC-tileplot.pdf'), width=5, height=4)
plot.GeTcD
dev.off()


plot.GeTcD <- ggplot(df.GeTcD, aes(GeT, TcD, z=auc_global)) +
  scale_fill_gradient(low='white', high='darkblue') +
  stat_contour(aes(fill=..level..), geom='polygon', bins=15) +
  geom_vline(xintercept=threshold.GeT, color='darkgreen', linetype='dashed') + 
  geom_hline(yintercept=threshold.DcT, color='darkgreen', linetype='dashed') + 
  theme_bw() + guides(fill=guide_colorbar(title='GeTcD AUC')) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

pdf(file.path(directory, 'AUC-contourplot.pdf'), width=5, height=4)
plot.GeTcD
dev.off()


# tissue localization
df.GeTcD <- subset(auc.df, metapath %in% c('GeTcD'))
plot.GeTcD <- ggplot(df.GeTcD, aes(GeT, TcD, fill=auc_logreg)) +
  geom_tile() + scale_fill_gradient(low='white', high='darkblue') +
  geom_vline(xintercept=threshold.GeT, color='darkgreen', linetype='dashed') + 
  geom_hline(yintercept=threshold.DcT, color='darkgreen', linetype='dashed') + 
  theme_bw() + guides(fill=guide_colorbar(title='GeTcD AUC')) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

pdf(file.path(directory, 'AUC-tileplot-ridge.pdf'), width=5, height=4)
plot.GeTcD
dev.off()

sorted.df.GeTcD <- df.GeTcD[order(df.GeTcD$auc_global, decreasing=TRUE), ]
write.table(sorted.df.GeTcD, file.path(directory, 'GeTcD-performance.txt'), sep='\t', quote=FALSE, row.names=FALSE)


#gridExtra::grid.arrange(plot.GeT, plot.GeTcD, ncol=1)


