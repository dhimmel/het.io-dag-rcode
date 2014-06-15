library(ggplot2)
library(grid)
library(plyr)
library(reshape2)

options(stringsAsFactors=FALSE)
options(width=Sys.getenv('COLUMNS'))



project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
prewtc.id <- '140522-no-wtccc2'
global.id <- '140522-all-assoc-lessmsig'

# Read Vegas
vegas.path <- file.path(project.dir, 'vegas', 'vegas_meta2.5_wtccc2_combined.txt')
vegas.df <- read.delim(vegas.path)


# Read Predictions
predictions.pre.wtc.path <- file.path(project.dir, 'networks', prewtc.id, 'MS-predictions-pre-wtccc2.txt')
predictions.pre.wtc.df <- read.delim(predictions.pre.wtc.path, check.names=FALSE)
predictions.pre.wtc.df <- predictions.pre.wtc.df[, c('gene_symbol', 'prediction')]
colnames(predictions.pre.wtc.df) <- c('gene', 'prior_pre_wtc')

predictions.post.wtc.path <- file.path(project.dir, 'networks', global.id, 'modeling', 'prediction-table.txt')
predictions.post.wtc.df <- read.delim(predictions.post.wtc.path, check.names=FALSE)
predictions.post.wtc.df <- predictions.post.wtc.df[, c('gene_symbol', 'multiple sclerosis')]
colnames(predictions.post.wtc.df) <- c('gene', 'prior_post_wtc')

prior.df <- merge(predictions.pre.wtc.df, predictions.post.wtc.df)

prior.scatplot <- ggplot(prior.df, aes(prior_pre_wtc, prior_post_wtc)) +
  geom_point(alpha=0.3) + geom_smooth() + 
  coord_trans(x='log10', y='log10') + theme_bw() +
  xlab('Pre-WTCCC2') + ylab('Post-WTCCC2') + ggtitle('MS Predictions') +
  theme(plot.margin = grid::unit(c(2, 2, 2, 2), 'points'))
path <- file.path(project.dir, 'MS-analysis', 'prior-scatterplot.pdf')
ggsave(path, prior.scatplot, width=4, height=6)

# 
ms.df <- merge(vegas.df, prior.df)

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

