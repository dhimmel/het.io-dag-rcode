library(ggplot2)
library(gridExtra)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
directory <- file.path(project.dir, 'networks', '140615-thresholdsweep')
aucs.path <- file.path(directory, 'feature-aucs.txt')
auc.df <- read.delim(aucs.path)


code.dir <- file.path(project.dir, 'rcode')
source(file.path(code.dir, 'functions.R'))



#threshold.DsD <- 0.2
#threshold.GfG <- 0.2
threshold.GeT <- 1.4
threshold.DlT <- 29

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
df.DlT <- subset(auc.df, metapath %in% c('GaDlTlD'))
plot.DlT <- ggplot(df.DlT, aes(DlT, auc_logreg, shape=metapath)) +
  geom_vline(xintercept=threshold.DlT, color='darkgreen', linetype='dashed') + 
  geom_point() + geom_line() + theme_bw() + ylab('AUC') +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

pdf(file.path(directory, 'AUC-lineplots-ridge.pdf'), width=4, height=4)
gridExtra::grid.arrange(plot.GeT, plot.DlT, ncol=1)
dev.off()


#### tissue localization


# Countourplot
plot.GeTlD <- ggplot(df.GeTlD, aes(GeT, TlD, z=auc_global)) +
  scale_fill_gradient(low='white', high='darkblue') +
  stat_contour(aes(fill=..level..), geom='polygon', bins=15) +
  geom_vline(xintercept=threshold.GeT, color='darkgreen', linetype='dashed') + 
  geom_hline(yintercept=threshold.DlT, color='darkgreen', linetype='dashed') + 
  theme_bw() + guides(fill=guide_colorbar(title='GeTlD AUC')) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

pdf(file.path(directory, 'AUC-contourplot.pdf'), width=5, height=4)
plot.GeTlD
dev.off()


# tileplot - ridge
df.GeTlD <- subset(auc.df, metapath %in% c('GeTlD'))
plot.GeTlD <- ggplot(df.GeTlD, aes(GeT, TlD, fill=auc_logreg)) +
  geom_tile() + scale_fill_gradient(low='white', high='darkblue') +
  geom_vline(xintercept=threshold.GeT, color='darkgreen', linetype='dashed') + 
  geom_hline(yintercept=threshold.DlT, color='darkgreen', linetype='dashed') + 
  theme_bw() + guides(fill=guide_colorbar(title='GeTlD AUC')) +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))

pdf(file.path(directory, 'AUC-tileplot-ridge.pdf'), width=5, height=4)
plot.GeTlD
dev.off()

sorted.df.GeTlD <- df.GeTlD[order(df.GeTlD$auc_global, decreasing=TRUE), ]
write.table(sorted.df.GeTlD, file.path(directory, 'GeTlD-performance.txt'), sep='\t', quote=FALSE, row.names=FALSE)


#gridExtra::grid.arrange(plot.GeT, plot.GeTlD, ncol=1)





# Thresholds - ridge
df.GeTlD <- subset(auc.df, metapath %in% c('GeTlD'))
plot.GeTlD <- ggplot(df.GeTlD, aes(GeT, TlD, fill=auc_global))
plot.GeTlD <- SetGGTheme(plot.GeTlD) +
  geom_tile() + 
  scale_fill_gradient(low='#F8F6F7', high='#300A24', breaks=seq(0.5, 1, 0.04)) +
  geom_vline(xintercept=threshold.GeT, color=Solar('yellow'), linetype='dashed') + 
  geom_hline(yintercept=threshold.DlT, color=Solar('yellow'), linetype='dashed') + 
  guides(fill=guide_colorbar(title='AUROC', barwidth=0.65, barheight=4, nbin=500)) +
  theme(legend.justification = c(1, 1), legend.position = c(1, 1)) +
  theme(legend.background=element_rect(color='grey60', size=0.2)) +
  scale_x_continuous(expand=c(0, 0)) + scale_y_continuous(expand=c(0, 0)) +
  xlab('Expression Threshold') + ylab('Localization Threshold')

pdf(file.path(directory, 'AUC-tileplot.pdf'), width=5, height=4)
plot.GeTlD
dev.off()








## Metrics
metric.path <- file.path(project.dir, 'networks', '140615-metricsweep', 'model-auc-10x20-CV.txt')
metric.df <- read.delim(metric.path, stringsAsFactors=FALSE)

metric.summary.df <- plyr::ddply(metric.df, c('alpha', 'metric', 'repitition'),
  plyr::summarize, 'auroc'=mean(auroc))
metric.summary.df[, 'dwpc_exponent'] <- suppressWarnings(as.numeric(gsub('DWPC_', '', metric.summary.df$metric)))

npc.ggdf <- subset(metric.summary.df, metric == 'NPC')
npc.ggdf <- plyr::ddply(npc.ggdf, 'alpha', plyr::summarize,
  'lower'=t.test(auroc, conf.level=0.9999)$conf.int[1],
  'upper'=t.test(auroc, conf.level=0.9999)$conf.int[2],
  'auroc'=mean(auroc))


threshold.dwpc <- 0.4

dwpc.ggdf <- subset(metric.summary.df, ! is.na(dwpc_exponent))
plot.alphas <- unique(dwpc.ggdf$alpha)
vline.df <- data.frame('threshold'=c(threshold.dwpc, rep(NA, length(plot.alphas) - 1)), 'alpha'=plot.alphas)

gg.metric <- ggplot(dwpc.ggdf, aes(dwpc_exponent, auroc)) 
gg.metric <- SetGGTheme(gg.metric) +
  facet_wrap(~ alpha, nrow=2) +
  geom_vline(data=vline.df, aes(xintercept=threshold),
    color=Solar('yellow'), linetype='dashed') + 
  geom_rect(data=npc.ggdf, aes(x=NULL, ymin=lower, ymax=upper), xmin=-Inf, xmax=Inf,
    fill=Solar('violet'), alpha=0.4) +
  geom_hline(data=npc.ggdf, aes(x=NULL, yintercept=auroc), color=Solar('violet'), size=0.25) +
  geom_smooth(linetype=0, method='loess', span=0.75, alpha=0.6,
    fill=Solar('magenta'), level=0.9999) +
  stat_summary(fun.y='mean', geom='point', color=Solar('magenta')) + 
  scale_x_continuous(breaks=seq(0.1, 0.7, 0.2)) +
  scale_y_continuous(breaks=seq(0.76, 0.80, 0.02)) + 
  #coord_cartesian(ylim=c(0.784, 0.821)) +
  theme(strip.text.x = element_text(size=9)) +
  xlab('DWPC Exponent') + ylab(expression(10 %*% 20 ~ fold ~ CV ~ AUROC))
#gg.metric

ExpressionFromAlpha <- function(a) {
  if (a == 0) {
    method.name <- '(ridge)'
  } else if (a == 1) {
    method.name <- '(lasso)'
  } else {
    method.name <- '(elastic net)'
  }
  expr.call <- bquote(alpha == .(a) ~ .(method.name))
  return(expr.call)
}

label.list <- lapply(unique(dwpc.ggdf$alpha), ExpressionFromAlpha)
gg.metric <- FacetWrapLabeller(gg.metric, label.list)


pdf(file.path(directory, 'optimization.pdf'), width=width.full, height=2.5)
gridExtra::grid.arrange(gg.metric, plot.GeTlD, nrow=1, widths=c(1.625, 1))
dev.off()






enet.auc.plot.path <- file.path(directory, sprintf('metrics-%sx%s-CV.pdf', repititions, kfolds))
ggsave(enet.auc.plot.path, gg.metric, width=width.half, height=3.6)

