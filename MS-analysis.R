options(stringsAsFactors=FALSE)


project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'

# Load code
code.dir <- file.path(project.dir, 'rcode')
sources <- c('settings.R', 'plotting.R')
for (source.filename in sources) {
  source(file.path(code.dir, source.filename))
}



###################################################################
### Histo

# Read Vegas
vegas.path <- '/home/dhimmels/Documents/serg/gwas/MS/vegas-cast.txt'
vegas.df <- read.delim(vegas.path, stringsAsFactors=FALSE)

hist.df <- cbind(vegas.df, strip_label='Meta2.5')
gg.hist <- ggplot(hist.df, aes(meta_p_value))
gg.hist <- SetGGTheme(gg.hist) +
  facet_grid(strip_label ~ .) +
  geom_hline(yintercept=1, color='grey') + 
  geom_histogram(aes(y = ..density..), breaks=seq(0, 1, 0.05), 
    fill=c(Solar('red'), Solar(rep('blue', 19))), color='white') +
  scale_x_continuous(breaks=breaks.prc, expand=c(0.03, 0)) + 
  scale_y_continuous(breaks=seq(0.2, 3, 0.4), expand=c(0.03, 0)) +
  xlab('P-value') + ylab('Density')

path <- file.path(project.dir, 'MS-analysis', 'meta2.5-vegas.pdf')
OpenPDF(path, width=width.half, height=2)
print(gg.hist)
ClosePDF(path)


#########################################################################

prewtc.id <- '140615-no-wtccc2'
global.id <- '140615-all-assoc'


# Read Predictions
predictions.pre.wtc.path <- file.path(project.dir, 'networks', prewtc.id, 'model', 'MS-predictions.txt.gz')
predictions.pre.wtc.df <- read.delim(predictions.pre.wtc.path, check.names=FALSE)
predictions.pre.wtc.df <- predictions.pre.wtc.df[, c('gene_symbol', 'ridge')]
colnames(predictions.pre.wtc.df) <- c('gene_symbol', 'prediction_pre_wtc')

predictions.post.wtc.path <- file.path(project.dir, 'networks', global.id, 'modeling', 'prediction-table.txt')
predictions.post.wtc.df <- read.delim(predictions.post.wtc.path, check.names=FALSE)
predictions.post.wtc.df <- predictions.post.wtc.df[, c('gene_symbol', 'multiple sclerosis')]
colnames(predictions.post.wtc.df) <- c('gene_symbol', 'prediction_post_wtc')

prior.df <- merge(predictions.pre.wtc.df, predictions.post.wtc.df)


# Prior Scatterplot
prior.scatplot <- ggplot(prior.df, aes(prediction_pre_wtc, prediction_post_wtc)) +
  geom_point(alpha=0.3) + geom_smooth() + 
  coord_trans(x='log10', y='log10') + theme_bw() +
  xlab('Pre-WTCCC2') + ylab('Post-WTCCC2') + ggtitle('MS Predictions') +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))
path <- file.path(project.dir, 'MS-analysis', 'prior-scatterplot.pdf')
ggsave(path, prior.scatplot, width=4, height=6)

# Merge Vegas and Predictions
vegas.select.df <- vegas.df[c('gene_symbol', 'chromosome', 'meta_base_start', 'meta_base_stop', 'meta_p_value', 'wtc_p_value')]
ms.df <- merge(vegas.select.df, prior.df)


# Read GWAS Catalog Associations
# Pre-WTCCC2
assoc.pre.wtc.path <- '/home/dhimmels/Documents/serg/data-sources/gwas-catalog/140205/processed-no-wtccc2/association-statuses.txt'
assoc.pre.wtc.df <- read.delim(assoc.pre.wtc.path, stringsAsFactors=FALSE)
assoc.pre.wtc.df <- subset(assoc.pre.wtc.df, disease_name == 'multiple sclerosis')
status.pre.wtc <- assoc.pre.wtc.df[match(ms.df$gene_symbol, assoc.pre.wtc.df$gene_symbol), 'status']
status.pre.wtc[is.na(status.pre.wtc)] <- 'negative'
ms.df$status_pre_wtc <- status.pre.wtc

# Post-WTCCC2
assoc.post.wtc.path <- '/home/dhimmels/Documents/serg/data-sources/gwas-catalog/140205/processed/association-statuses.txt'
assoc.post.wtc.df <- read.delim(assoc.post.wtc.path, stringsAsFactors=FALSE)
assoc.post.wtc.df <- subset(assoc.post.wtc.df, disease_name == 'multiple sclerosis')
status.post.wtc <- assoc.post.wtc.df[match(ms.df$gene_symbol, assoc.post.wtc.df$gene_symbol), 'status']
status.post.wtc[is.na(status.post.wtc)] <- 'negative'
ms.df$status_post_wtc <- status.post.wtc

# status_pre_wtc status_post_wtc
hcp_post_wtc <- subset(ms.df, status_post_wtc == 'HC_primary')[, 'gene_symbol']
hcs_post_wtc <- subset(ms.df, status_post_wtc == 'HC_secondary')[, 'gene_symbol']

ms.df <- ms.df[order(ms.df$chromosome, ms.df$meta_base_start, ms.df$meta_base_stop), ]
path <- file.path(project.dir, 'MS-analysis', 'MS-predictions.txt')
write.table(ms.df, path, quote=FALSE, row.names=FALSE, sep='\t')

# xMHC region defined from doi:10.1038/nrg1489
xmhc <- ms.df[which(ms.df$gene_symbol == 'SCGN'):which(ms.df$gene_symbol == 'SYNGAP1'), 'gene_symbol']

# WTCCC2 Replicated genes http://wattle.well.ox.ac.uk/wtccc2/external/ms/
path <- '/home/dhimmels/Documents/serg/gwas/MS/WTCCC2/genes-that-replicated.txt'
wtc.repl.df <- read.table(path, sep='\t', col.names=c('primary', 'secondary'), fill=TRUE, na.strings='')
wtc.replicated <- as.character(na.omit(unlist(wtc.repl.df)))
wtc.replicated <- wtc.replicated[wtc.replicated %in% ms.df$gene_symbol]

# Calculated linked genes
linked <- unique(c(wtc.replicated, hcp_post_wtc))
for (primary in linked) {
  i <- which(ms.df$gene_symbol == primary)
  chrom.df <- subset(ms.df, chromosome == ms.df[i, 'chromosome'])
  i <- which(chrom.df$gene_symbol == primary)
  is.nominal <- chrom.df$wtc_p_value <= 0.05
  upstream <- cumsum(is.nominal[i:1]) == seq.int(i)
  downstream <- cumsum(is.nominal[i:nrow(chrom.df)]) == seq.int(nrow(chrom.df) - i + 1)
  linked_symbols <- chrom.df[c(upstream, downstream[-1]), 'gene_symbol']
  linked <- c(linked, linked_symbols)
}
exclusions <- unique(c(linked, hcs_post_wtc, xmhc))
novel.df <- subset(ms.df, ! (gene_symbol %in% exclusions))
novel.df <- novel.df[order(novel.df$prediction_pre_wtc, decreasing=TRUE), ]
nominal.novel.df <- subset(novel.df, meta_p_value <= 0.05) 

BonferroniValidator <- function(pvals, alpha=0.05) {
  bv.rows <- lapply(1:length(pvals), function(i) {
    cutoff <- alpha / i
    q <- sum(pvals[1:i] < cutoff)
    m <- sum(pvals < cutoff)
    n <- length(pvals) - m
    pval <- phyper(q - 1, m, n, i, lower.tail=FALSE)
    return(c(cutoff, q, pval))
  })
  bv.df <- as.data.frame(do.call(rbind, bv.rows))
  colnames(bv.df) <- c('bonferroni_cutoff', 'number_validated', 'p_value_hyper')
  return(bv.df)
}

bv.df <- BonferroniValidator(nominal.novel.df$wtc_p_value)
nominal.novel.df <- cbind(nominal.novel.df, bv.df)
path <- file.path(project.dir, 'MS-analysis', 'MS-novel-predictions.txt')
write.table(nominal.novel.df, path, quote=FALSE, row.names=FALSE, sep='\t')

