
options(stringsAsFactors=FALSE)
options(width=Sys.getenv('COLUMNS'))

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'

# Read Meta2.5 Vegas
meta.vegas.path <- file.path(project.dir, 'vegas', 'meta2.5-genewise.txt')
meta.vegas.df <- read.delim(meta.vegas.path, na.strings=c('#N/A', ''))
meta.vegas.df <- data.frame(
  'gene'=meta.vegas.df$Gene,
  'chromosome'=meta.vegas.df$Chr,
  'location_start'=meta.vegas.df$Start,
  'location_stop'=meta.vegas.df$Stop,
  'pval_meta'=meta.vegas.df$Pvalue_all,
  'mlogp_meta'=-log10(meta.vegas.df$Pvalue_all)
)


# Read WTCCC2 Vegas
wtc.vegas.path <- file.path(project.dir, 'vegas', 'wtccc2-genewise.txt')
wtc.vegas.df <- read.delim(wtc.vegas.path, na.strings=c('#N/A', ''))
wtc.genes <- wtc.vegas.df[as.logical(wtc.vegas.df[, 'Is.Nature.gene']), 'Gene']
wtc.blocks <- wtc.vegas.df[as.logical(wtc.vegas.df[, 'Is.Nature.gene']), 'block']
manual.blocks <- c(226, 227, 230, 232, 233) # should we change to entire region
wtc.blocks <- append(wtc.blocks, manual.blocks)
wtc.blocks <- unique(na.omit(wtc.blocks))
wtc.linked.genes <- wtc.vegas.df[wtc.vegas.df$block %in% wtc.blocks, 'Gene']
wtc.linked.genes <- unique(c(wtc.genes, wtc.linked.genes))
wtc.vegas.df$wtc.linked <- as.integer(wtc.vegas.df$Gene %in% wtc.linked.genes)
wtc.vegas.df <- data.frame(
  'gene'=wtc.vegas.df$Gene,
  'pval_wtc'=wtc.vegas.df$P.value_all.SNPs,
  'mlogp_wtc'=-log10(wtc.vegas.df$P.value_all.SNPs),
  'wtc_gene'=wtc.vegas.df$Is.Nature.gene,
  'wtc_linked'=wtc.vegas.df$wtc.linked
)


# Merge
vegas.df <- merge(meta.vegas.df, wtc.vegas.df)
vegas.df <- vegas.df[order(vegas.df$chromosome, vegas.df$location_start), ]
vegas.df <- subset(vegas.df, ! (is.na(pval_meta) | is.na(pval_wtc)))

vegas.merged.path <-file.path(project.dir, 'vegas', 'vegas_meta2.5_wtccc2_combined.txt')
write.table(vegas.df, vegas.merged.path, sep='\t', row.names=FALSE, quote=FALSE, na='')

