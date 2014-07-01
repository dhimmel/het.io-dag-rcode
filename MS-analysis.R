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

linked <- NULL
for (primary in hcp_post_wtc) {
  i <- which(ms.df$gene_symbol == primary)
  chrom.df <- subset(ms.df, chromosome == ms.df[i, 'chromosome'])
  i <- which(chrom.df$gene_symbol == primary)
  is.nominal <- chrom.df$wtc_p_value <= 0.05
  upstream <- cumsum(is.nominal[i:1]) == seq.int(i)
  downstream <- cumsum(is.nominal[i:nrow(chrom.df)]) == seq.int(nrow(chrom.df) - i + 1)
  linked_symbols <- chrom.df[c(upstream, downstream[-1]), 'gene_symbol']
  linked <- c(linked, linked_symbols)
}
exclusions <- unique(c(hcp_post_wtc, linked, hcs_post_wtc))
novel.df <- subset(ms.df, ! (gene_symbol %in% exclusions))
novel.df <- novel.df[order(novel.df$prediction_pre_wtc, decreasing=TRUE), ]
nominal.novel.df <- subset(novel.df, meta_p_value <= 0.05) 



novel.df <- subset(ms.df, ! wtc_linked)
nominal.novel.df <- subset(novel.df, pval_meta <= 0.05)
nominal.novel.df <- nominal.novel.df[order(nominal.novel.df$prior_pre_wtc, decreasing=TRUE), ]


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
  colnames(bv.df) <- c('bonf.cutoff', 'num.validated', 'pval.hypge')
  return(bv.df)
}

bv.df <- BonferroniValidator(nominal.novel.df$pval_wtc)
nominal.novel.df <- cbind(nominal.novel.df, bv.df)

nominal.novel.df$validated_genes <- sapply(1:nrow(nominal.novel.df), function(i) {
  cutoff <- nominal.novel.df[i, 'bonf.cutoff']
  subset.df <- nominal.novel.df[1:i, ]
  gene.vec <- subset.df[subset.df$pval_wtc < cutoff, 'gene']
  paste(gene.vec, collapse='|')
})

path <- file.path(project.dir, 'MS-analysis', 'nominal-in-meta2.5-novel.txt')
write.table(nominal.novel.df, path, sep='\t', row.names=FALSE, quote=FALSE)



















#dapple.df <- read.delim(file.path(project.dir, 'comparison', 'dapple', 'MSnowtccc2genes2_CIscores.txt'))
#dapple.df$mlog_pval <- -log10(dapple.df$P_VALUE)
#dapple.df <- subset(dapple.df, PROTEIN %in% features.df$source)
#features.df$dapple <- 0
#features.df[match(dapple.df$PROTEIN, features.df$source), 'dapple'] <- dapple.df$mlog_pval

# Read immunochip genes
immunochip.df <- read.delim(file.path(project.dir, 'data-integration', 'MS-immunochip-novel-variants.txt'), na.strings='')
genes.immunochip.novel <- unique(na.omit(immunochip.df$gene))
predictions.wtc.df <- read.delim('/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/140321-all-assoc/prediction-table.txt', check.names=FALSE)
ms.df$prediction_with_wtc <- predictions.wtc.df[match(ms.df$source, predictions.wtc.df$source), 'multiple sclerosis']
ms.df$ichip_novel <- as.integer(ms.df$source %in% genes.immunochip.novel)

ms.novel.ichip.df <- subset(ms.df, (! status) & (! wtccc2_novel))

vtm.ichip <- VariableThresholdMetrics(ms.novel.ichip.df$prediction_with_wtc, ms.novel.ichip.df$ichip_novel)



#################################################
## VEGAS


wtc.df <- read.delim(file.path(project.dir, 'vegas', 'wtccc2-genewise.txt'), na.strings=c('#N/A', ''))
meta.df <- meta.df[, c('Gene', 'Pvalue_all')]
wtc.df <- wtc.df[, c('Gene', 'P.value_all.SNPs')]
colnames(meta.df) <- c('source', 'pval.meta')
colnames(wtc.df) <- c('source', 'pval.wtc')
vegas.df <- merge(meta.df, wtc.df)
vegas.df$mlogp.meta <- -log10(vegas.df$pval.meta)
vegas.df$mlogp.wtc <- -log10(vegas.df$pval.wtc)
vegas.df <- merge(vegas.df, ms.df[, c('source', 'prediction', 'status', 'wtccc2_novel')])

RankValidator <- function(score, validation.p) {
  ranked.p <- validation.p[order(score, decreasing=TRUE)]
  sapply(1:length(ranked.p), function(i) sum(ranked.p[1:i] * i < 0.05, na.rm=TRUE))
}

CalculatePPA <- function(prior, pval) {
  # http://www.nature.com/nrg/journal/v10/n10/full/nrg2615.html
  # BF - Bayes Factor
  # PO - Posterior Odds
  # PPA - Posterior Probability of Association
  BF <- -1 / (exp(1) * pval * log(pval))
  BF[pval > exp(-1)] <- NA
  PO <- BF * prior / (1 - prior)
  PPA <- PO / (1 + PO)
  return(PPA)
}

vegas.df$posterior <- CalculatePPA(vegas.df$prediction, vegas.df$pval.meta)

vegas.df <- subset(vegas.df, pval.meta < 0.05 & pval.meta > 10e-5)

rank.df <- data.frame(
  'random'=RankValidator(sample(1:nrow(vegas.df)), vegas.df$pval.wtc),
  'prior'=RankValidator(vegas.df$prediction, vegas.df$pval.wtc),
  'meta2.5'=RankValidator(vegas.df$mlogp.meta, vegas.df$pval.wtc),
  'posterior'=RankValidator(vegas.df$posterior, vegas.df$pval.wtc))
rank.df$rank <- 1:nrow(rank.df)
rank.melt <- reshape2::melt(rank.df, id.vars=c('rank'), variable.name='method', value.name='validated')


ggplot(rank.melt, aes(rank, validated, color=method)) +
  geom_point() + theme_bw()


#plot(vegas.df$mlogp.meta, vegas.df$mlogp.wtc)
#ggplot(vegas.df, aes(mlogp.meta, mlogp.wtc)) + geom_point() + geom_smooth()

