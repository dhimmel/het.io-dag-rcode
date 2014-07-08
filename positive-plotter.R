options(stringsAsFactors=FALSE)

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
code.dir <- file.path(project.dir, 'rcode')
sources <- c('settings.R', 'plotting.R')
for (source.filename in sources) {
  source(file.path(code.dir, source.filename))
}

path <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/140615-training/model/testing-predictions-ridge.txt.gz'
ridge.df <- read.delim(path, stringsAsFactors=FALSE, check.names=FALSE)
ridge.df$prediction <- 100 * ridge.df$prediction

path <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/permuted/140615-0-training/model/testing-predictions-ridge.txt.gz'
perm.df <- read.delim(path, stringsAsFactors=FALSE, check.names=FALSE)
perm.df$prediction <- 100 * perm.df$prediction


PercentDF <- function(prediction.df) {
  pred.neg <- subset(prediction.df, status == 0)$prediction
  pred.pos <- subset(prediction.df, status == 1)$prediction
  quantiles <- quantile(pred.pos, seq(0, 1, 0.1), dig.lab=10)
  cuts.pos <- cut(pred.pos, breaks=quantiles)
  cuts.neg <- cut(pred.neg, breaks=quantiles)
  table.neg <- table(cuts.neg)
  table.pos <- table(cuts.pos)
  percent.df <- data.frame('quantile'=levels(cuts.pos), 
    'positives'=as.numeric(table.pos), 
    'negatives'=as.numeric(table.neg))
  percent.df$total <- percent.df$positives + percent.df$negatives
  percent.df$percent_positive <- percent.df$positives / percent.df$total
  return(percent.df)
}


percent.df <- rbind(
  cbind(PercentDF(ridge.df), 'panel'='Ridge'),
  cbind(PercentDF(perm.df), 'panel'='Permuted Ridge'))

ratio.df <- subset(percent.df, panel == 'Ridge')
ratio <- subset(percent.df, panel == 'Ridge')$percent_positive /
  subset(percent.df, panel == 'Permuted Ridge')$percent_positive
ratio.df$ratio <- ratio
ratio.df$fratio <- format(ratio, digits=2)

gg.bar <- ggplot(percent.df, aes(x=quantile, y=percent_positive))
gg.bar <- SetGGTheme(gg.bar) +
  facet_grid(. ~ panel, scales='free_x') +
  geom_bar(stat='identity', fill=Solar('base01')) +
  theme(axis.text.x=element_text(angle=35, hjust=1)) +
  xlab('Prediction (%) Quantile') + ylab('Percent Positive') +
  geom_text(data=ratio.df, aes(label=fratio), size=3.1, vjust=-0.2) +
  ylim(c(0, max(percent.df$percent_positive * 1.065)))
gg.bar

path <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/permuted/140615-0-training/plots/prediction-quantiles.pdf'
OpenPDF(path, width=width.full, height=3)
print(gg.bar)
ClosePDF(path)

path <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/permuted/140615-0-training/model/prediction-quantiles.txt'
write.table(percent.df, path, row.names=FALSE, quote=FALSE, sep='\t')
path <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/permuted/140615-0-training/model/prediction-quantile-ratios.txt'
write.table(ratio.df, path, row.names=FALSE, quote=FALSE, sep='\t')

