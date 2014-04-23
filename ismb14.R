
ismb.df <- rbind(
data.frame(
'disease'=disease.auc.df$target_name, 
'category'=disease.auc.df$category,
'auc'=disease.auc.df$auc,
'xname'='Global',
'type'='Model'),
data.frame(
'disease'=msig.auc.df$target_name, 
'category'=msig.auc.df$category,
'auc'=msig.auc.df$auc,
'xname'=msig.auc.df$short_name,
'type'='Gene-{MSigDB Collection}-Gene-Disease Feature')
)

global.ismb.df <- rbind(
data.frame('type'='Model', 'xname'='Global', 'auc'=vtm.global$auc),
data.frame(
'type'='Gene-{MSigDB Collection}-Gene-Disease Feature',
'xname'=global.msig.auc.df$short_name,
'auc'=global.msig.auc.df$auc)
)

global.msig.auc.df$short_name <- as.character(global.msig.auc.df$short_name)
global.msig.auc.df$short_name[global.msig.auc.df$short_name == 'Cancer Neighborhoods'] <- 'Cancer Hoods'
global.ismb.df$xname[global.ismb.df$xname == 'Cancer Neighborhoods'] <- 'Cancer Hoods'
ismb.df$xname[ismb.df$xname == 'Cancer Neighborhoods'] <- 'Cancer Hoods'

msig.levels <- c(as.character(global.msig.auc.df$short_name[order(global.msig.auc.df$auc)]), 'Global')

ismb.df$xname <- factor(ismb.df$xname, levels=msig.levels)
global.ismb.df$xname <- factor(global.ismb.df$xname, levels=msig.levels)


ismb.plot <- ggplot(ismb.df, aes(xname, auc)) + 
  facet_grid(. ~ type, scales='free_x', space='free_x') +
  geom_hline(aes(yintercept=0.5), color='red', linetype='dashed') +
  stat_summary(fun.y='mean.conf.int', geom='line', color='grey', size=12) +
  geom_point(data=global.ismb.df, size=15, shape='-', color='#4C4C4C') + 
  geom_point(aes(color=category), position=position_jitter(width=0.3), alpha=0.7, size=2.5) + 
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(legend.position='top') +
  xlab(NULL) + ylab('AUROC') + 
  scale_colour_manual(name='Disease Category', values=category.colors)
path <- file.path(network.dir, 'ismb-plot.pdf')
ggsave(path, ismb.plot, width=9, height=5)































##############################################################################
#### Old AUROC Plotting





### WE ARE HERE

category.path <- file.path(project.dir, 'data-integration', 'disease-categories.txt')
category.df <- read.delim(category.path)
category.merge.df <- data.frame('target'=category.df$doid_code, 'category'=category.df$category)
feature.auc.melt <- merge(feature.auc.melt, category.merge.df)
feature.auc.melt <- subset(feature.auc.melt, feature != 'PC_t|G-a-D')




feature.auc.melt$metapath.abbrev <- gsub('-', '', gsub('DWPC_0.4\\|', '', feature.auc.melt$feature))
global.feature.auc.df$metapath.abbrev <- gsub('-', '', gsub('DWPC_0.4\\|', '', global.feature.auc.df$feature))


path <- file.path(network.dir, 'feature-aucs-by-disease.txt')
write.table(feature.auc.melt, path, sep='\t', row.names=FALSE, quote=FALSE)
path <- file.path(network.dir, 'feature-aucs-global.txt')
write.table(global.feature.auc.df, path, sep='\t', row.names=FALSE, quote=FALSE)




##
# MsigDB Graphic
msig.nomen.path <- file.path(project.dir, 'data-integration', 'MSigDB-nomenclature.txt')
msig.nomen.df <- read.delim(msig.nomen.path)
msig.nomen.df$short_name <- factor(msig.nomen.df$short_name, levels=msig.nomen.df$short_name)
msig.nomen.df$metapath.abbrev <- paste('Gm', msig.nomen.df$code, 'mGaD', sep='')

msig.auc.df <- merge(feature.auc.melt, msig.nomen.df[, c('metapath.abbrev', 'short_name')])
global.msig.auc.df <- merge(global.feature.auc.df, msig.nomen.df[, c('metapath.abbrev', 'short_name')])

auc.msig.plot <- ggplot(msig.auc.df, aes(short_name, auc)) + 
  geom_hline(aes(yintercept=0.5), color='red', linetype='dashed') +
  stat_summary(fun.y='mean.conf.int', geom='line', color='grey', size=10) +
  geom_point(data=global.msig.auc.df, size=13, shape='-', color='#4C4C4C') + 
  geom_point(aes(color=category), position=position_jitter(width=0.2), alpha=0.7, size=2.5) + 
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + theme(legend.position='top') +
  xlab('DWPC for the Gene-{Gene Set}-Gene-Disease Metapath') + ylab('AUROC') +
  scale_colour_manual(name='Disease Category', values=category.colors)
path <- file.path(network.dir, 'AUC-msig-plot.pdf')
ggsave(path, auc.msig.plot, width=7.5, height=5)

##
# Non-MsigDB Graphic

no.msig.feature.auc.df <- subset(feature.auc.melt, ! metapath.abbrev %in% msig.nomen.df$metapath.abbrev)
no.msig.global.feature.auc.df <- subset(global.feature.auc.df, ! metapath.abbrev %in% msig.nomen.df$metapath.abbrev)

ggplot(no.msig.feature.auc.df, aes(metapath.abbrev, auc)) + 
  geom_hline(aes(yintercept=0.5), color='red', linetype='dashed') +
  stat_summary(fun.y='mean.conf.int', geom='line', color='grey', size=10) +
  geom_point(data=no.msig.global.feature.auc.df, size=20, shape='-', color='#4C4C4C') + 
  geom_point(aes(color=category), position=position_jitter(width=0.2), alpha=0.5) + 
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  xlab('Metapath') + ylab('AUROC') + theme(legend.position='top') +
  scale_colour_manual(name='Disease Category', values=category.colors)




# Full Model Graphic
disease.auc.df <- merge(disease.auc.df, category.merge.df)
disease.auc.df$model <- 'Full Model'
ggplot(disease.auc.df, aes(model, auc)) + 
  geom_hline(aes(yintercept=0.5), color='red', linetype='dashed') +
  stat_summary(fun.y='mean.conf.int', geom='line', color='grey', size=10) +
  geom_point(data=data.frame('model'='Full Model', 'auc'=vtm.global$auc), size=20, shape='-', color='#4C4C4C') + 
  geom_point(aes(color=category), position=position_jitter(width=0.2), alpha=0.5) + 
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  xlab('Model') + ylab('AUROC') + 
  scale_colour_manual(name='Disease Category', values=category.colors)


no.msig.feature.auc.df$model <- no.msig.feature.auc.df$metapath.abbrev
no.msig.feature.auc.df$type <- 'DWPC Features'
no.msig.feature.auc.df[substr(as.character(no.msig.feature.auc.df$feature), 1, 2) == 'PC', 'type'] <- 'PC Features'
disease.auc.df$type <- 'Models'
no.msig.auc.df <- rbind(disease.auc.df[, c('target_name', 'category', 'model', 'type', 'auc')],
  no.msig.feature.auc.df[, c('target_name', 'category', 'model', 'type', 'auc')])


no.msig.global.feature.auc.df$type <- 'DWPC Features'
no.msig.global.feature.auc.df[substr(as.character(no.msig.global.feature.auc.df$feature), 1, 2) == 'PC', 'type'] <- 'PC Features'
no.msig.model.auc.df <- data.frame('model'=c('Global', no.msig.global.feature.auc.df$metapath.abbrev),
  'auc'=c(vtm.global$auc, no.msig.global.feature.auc.df$auc),
  'type'=c('Models', no.msig.global.feature.auc.df$type))


no.msig.auc.df$type <- factor(no.msig.auc.df$type, levels=c('DWPC Features', 'PC Features', 'Models'))
no.msig.model.auc.df$type <- factor(no.msig.model.auc.df$type, levels=c('DWPC Features', 'PC Features', 'Models'))

auc.plot <- ggplot(no.msig.auc.df, aes(model, auc)) + 
  facet_grid(. ~ type, scales='free_x', space='free_x') +
  geom_hline(aes(yintercept=0.5), color='red', linetype='dashed') +
  stat_summary(fun.y='mean.conf.int', geom='line', color='grey', size=12) +
  geom_point(data=no.msig.model.auc.df, size=15, shape='-', color='#4C4C4C') + 
  geom_point(aes(color=category), position=position_jitter(width=0.3), alpha=0.7, size=2.5) + 
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points')) + theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  theme(legend.position='top') +
  xlab(NULL) + ylab('AUROC') + 
  scale_colour_manual(name='Disease Category', values=category.colors)
path <- file.path(network.dir, 'AUC-plot.pdf')
ggsave(path, auc.plot, width=7.5, height=5)


size.plot <- ggplot(disease.auc.df, aes(positives, auc)) +
  geom_point(aes(color=category), size=4, alpha=0.7) + #geom_smooth() +
  scale_colour_manual(name='Disease Category', values=category.colors) +
  theme_bw() + xlab('Number of Positives') + ylab('AUROC')
path <- file.path(network.dir, 'AUC-versus-positives.pdf')
ggsave(path, size.plot, width=5, height=3)


