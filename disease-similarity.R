library(glmnet)
library(ggplot2)
library(doMC)
library(grid)
library(plyr)
library(reshape2)

options(stringsAsFactors=FALSE)
options(width=Sys.getenv('COLUMNS'))

project.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/'
network.dir <- file.path(project.dir, 'networks', '140615-all-assoc')

pred.tab.path <- file.path(network.dir, 'model', 'prediction-table-select.txt')
pred.tab <- read.delim(pred.tab.path, check.names=FALSE, row.names=1)

#pred.table <- log(pred.tab)
plot(density(pred.table[,2]))

cor.mat <- cor(pred.tab, method='spearman')
cor.mat[lower.tri(cor.mat, diag=TRUE)] <- NA
spearman.df <- na.omit(reshape2::melt(cor.mat))

#cor.path <- file.path(network.dir, 'spearman-disease-correlations.txt')
#write.table(spearman.df, cor.path, sep='\t', row.names=FALSE, quote=FALSE)

distance.mat <- 1 - cor(pred.tab, method='spearman')

#plot(hclust(as.dist(distance.mat), method='complete'))

plot(hclust(as.dist(distance.mat), method='ward'))


kmc <- kmeans(distance.mat, centers=2)


kmc <- kmeans(distance.mat, centers=5)
cbind(sort(kmc$cluster))




#image(cor.mat)


#http://weitaiyun.blogspot.com/2009/03/visulization-of-correlation-matrix.html
