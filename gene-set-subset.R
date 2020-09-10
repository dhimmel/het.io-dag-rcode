options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'

# Load code
code.dir <- file.path(project.dir, 'rcode')
sources <- c('settings.R', 'hnlp-learning.R', 'plotting.R')
for (source.filename in sources) {
  source(file.path(code.dir, source.filename))
}

library(dplyr)

# Load network info
network.id <- '140615-all-assoc'
network.dir <- file.path(project.dir, 'networks', network.id)
dirs <- InitializeNetworkDir(network.dir)

# Read the gene set subsetting performance results
auroc.df <- file.path(dirs$model, 'gene-set-subset-aurocs.txt') %>%
  read.delim()

min_nodes <- min(auroc.df$total_nodes)
min_edges <- min(auroc.df$total_edges)

auroc.df <- auroc.df %>%
  dplyr::mutate(x = ifelse(mask_type=='node', 100 * nodes / total_nodes, 100 * edges / total_edges)) %>%
  dplyr::mutate(subset_kind = ifelse(
    (mask_type=='node' & nodes == min_nodes) | (mask_type=='edge' & edges == min_edges),
    'minimum', 'range'))

# Create human readable node names
msigdb.df <- file.path(project.dir, 'data-integration', 'MSigDB-nomenclature.txt') %>%
  read.delim() %>%
  dplyr::select(metanode=filename_abbrev, metanode_name=short_name) %>%
  dplyr::bind_rows(data.frame(metanode='tissue', metanode_name='Tissue'))

msigdb.df <- auroc.df %>%
  dplyr::filter(x == 100) %>% 
  dplyr::group_by(metanode) %>% 
  dplyr::summarize(auroc=mean(auroc)) %>%
  dplyr::right_join(msigdb.df)  %>%
  dplyr::filter(! is.na(auroc))

msigdb.df$metanode_name <- factor(msigdb.df$metanode_name,
  levels=msigdb.df$metanode_name[order(msigdb.df$auroc, decreasing=TRUE)])

# label.df
label.df <- auroc.df %>%
  dplyr::mutate(total_nodes = ifelse(metanode == 'tissue', 77, total_nodes)) %>%
  dplyr::group_by(metanode) %>% 
  dplyr::summarize(
    mean_degree = sprintf('d: %s', format(sum(total_edges) / sum(total_nodes), big.mark=',', digits=3) ),
    total_nodes = sprintf('n: %s', format(mean(total_nodes), big.mark=',')), 
    total_edges = sprintf('e: %s', format(mean(total_edges), big.mark=','))) %>% 
  dplyr::left_join(msigdb.df)


lower.label.df <- label.df %>%
  dplyr::filter(metanode_name %in% levels(metanode_name)[1:5]) %>%
  dplyr::mutate(x=52, auroc=-Inf)

upper.label.df <- label.df %>%
  dplyr::filter(metanode_name %in% levels(metanode_name)[6:15]) %>%
  dplyr::mutate(x=0, auroc=Inf)

# plot
color.values <- c('node'=Solar('violet'), 'edge'=Solar('magenta'))
auroc.df$mask_type <- factor(auroc.df$mask_type, levels=c('node', 'edge'))

gg.auroc <- auroc.df %>%
  dplyr::filter(! (mask_type == 'edge' & subset_kind == 'minimum')) %>%
  dplyr::filter(! (metanode == 'tissue' & subset_kind == 'minimum')) %>%
  dplyr::left_join(msigdb.df %>% dplyr::select(metanode, metanode_name)) %>%
  ggplot(aes(x, auroc)) %>% SetGGTheme() +
  facet_wrap( ~ metanode_name, ncol = 5) +
  geom_smooth(aes(fill=mask_type), linetype=0, method = 'loess', alpha=0.6) +
  geom_point(aes(shape=subset_kind, alpha = subset_kind, color=mask_type)) +
  scale_shape_manual(values=c('range'=16, 'minimum'=4), guide=FALSE) +
  scale_alpha_manual(values=c('range'=0.2, 'minimum'=0.8), guide=FALSE) +
  scale_color_manual(values=color.values, name='Masked by') +
  scale_fill_manual(values=color.values, name='Masked by') +
  theme(legend.justification=c(1, 1), legend.position=c(1, 1),
    legend.background=element_rect(color = 'grey', size = 0.3)) +
  xlab('Percent of Nodes/Edges Masked') + ylab('AUROC') +
  scale_x_continuous(breaks=seq(10, 90, 20), labels=seq(90, 10, -20)) +

  # add stats to upper left
  geom_text(aes(label=total_nodes), data=lower.label.df, hjust=0, vjust=-2.9, size=3.3) +
  geom_text(aes(label=total_edges), data=lower.label.df, hjust=0, vjust=-1.7, size=3.3) +
  geom_text(aes(label=mean_degree), data=lower.label.df, hjust=0, vjust=-0.5, size=3.3) +
  # add stats to lower right
  geom_text(aes(label=total_nodes), data=upper.label.df, hjust=0, vjust=1.3, size=3.3) +
  geom_text(aes(label=total_edges), data=upper.label.df, hjust=0, vjust=2.5, size=3.3) +
  geom_text(aes(label=mean_degree), data=upper.label.df, hjust=0, vjust=3.7, size=3.3)

path <- file.path(dirs$plots, 'gene-set-subset-aurocs.pdf')
OpenPDF(path, width=7.5, height=8.75)
print(gg.auroc)
ClosePDF(path)

