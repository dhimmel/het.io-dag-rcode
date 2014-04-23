library(ggplot2)
library(gridExtra)

options(width=Sys.getenv('COLUMNS'))

assoc.df <- read.delim('/home/dhimmels/Documents/serg/data-sources/gwas-catalog/140205/associations.txt')

associations.per.gene <- table(assoc.df$symbol)
gene.df <- data.frame('symbol'=names(associations.per.gene), 'associations'=as.numeric(associations.per.gene))

ggplot(gene.df, aes(associations)) +
  geom_histogram(binwidth=1, alpha=0.5) + scale_y_log10() + theme_bw()


dbinom(1:6, size=nrow(assoc.df), prob=1/18000)
