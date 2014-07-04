library(ggplot2)
library(grid)
library(gtools)
library(gridExtra)

# Solarized Color Pallete
solarized <- c(
  yellow='#b58900', orange='#cb4b16', red='#dc322f', magenta='#d33682',
  violet='#6c71c4', blue='#268bd2', cyan='#2aa198', green='#859900',
  base03='#002b36', base02='#073642', base01='#586e75', base00='#657b83',
  base0='#839496', base1='#93a1a1', base2='#eee8d5', base3='#fdf6e3')

Solar <- function(...) {
  color.names <- c(...)
  hex.cols <- as.character(solarized[color.names])
  return(hex.cols)
}

# plotting parameters
width.full <- 6.83
width.half <- 3.27

breaks.roc <- seq(0, 1, by=0.2)
breaks.prc <- seq(0, 1, by=0.2)
ylim.prc <- c(0, 1)

# Colors
strip.fill <- Solar('base3')


'#300A24' # Blackberry (ubuntu terminal)
'#F8F6F7' # lighter blackberry



## Plotting functions

OpenPDF <- function(path, width=4, height=4) {
  cairo_pdf(path, width=width, height=height, bg='white', family='sans')
}

ClosePDF <- function(path) {
  dev.off()
  embedFonts(path)
}

ChrRound <- function(x, digits=2) {
  # Round x to digits and return as a character with trailing zeros.
  x.rounded <- round(x, digits)
  sprintf.chr <- sprintf('%%.%sf', digits)
  x.chr <- sprintf(sprintf.chr, x.rounded)
  return(x.chr)
}


SetGGTheme <- function(gg) {
  gg <- gg + theme_bw() + theme(legend.title=element_text(face='plain'))
  gg <- gg + theme(plot.margin=grid::unit(c(2, 2, 2, 2), 'points'))
  gg <- gg + theme(legend.margin=grid::unit(0, 'cm'))
  gg <- gg + theme(legend.key.height=grid::unit(0.90, 'lines'))
  gg <- gg + theme(strip.background=element_rect(fill=strip.fill))
  return(gg)
}


ggROC <- function(gg) {
  gg <- SetGGTheme(gg) + 
    geom_segment(x=0, xend=1, y=0, yend=1, color='grey', size=0.15, show_guide=FALSE) +
    theme(legend.background=element_rect(color='grey60', size=0.2)) +
    theme(legend.key.width=grid::unit(1.5, 'lines')) +
    theme(legend.key=element_rect(linetype='blank')) +
    theme(legend.justification=c(1, 0), legend.position=c(1, 0)) +
    scale_x_continuous(breaks=breaks.roc, expand=c(0.03, 0)) + 
    scale_y_continuous(breaks=breaks.roc, expand=c(0.03, 0)) +
    xlab('False Positive Rate') + ylab('Recall') + coord_fixed()
  return(gg)
}

ggPRC <- function(gg) {
  gg <- SetGGTheme(gg) +
  geom_path(size=.8, color='grey') + geom_point(size=1.5) +
  theme(legend.justification=c(1, 1), legend.position=c(1, 1)) +
  theme(legend.background=element_rect(color='grey60', size=0.2)) +
  xlab('Recall') + ylab('Precision') +
  scale_x_continuous(breaks=breaks.prc, expand=c(0.03, 0)) + 
  scale_y_continuous(breaks=breaks.prc, expand=c(0.03, 0), limits=ylim.prc) +
  #scale_color_gradientn(colours = rainbow(7)[6:1]) +
  scale_color_gradientn(colours = as.character(solarized[8:2])) +
  guides(color=guide_colorbar(title='Prediction\nThreshold', nbin=500))
  return(gg)
}

FacetWrapLabeller <- function(gg.plot, label.list=NULL) {
  #works with R 3.0.1 and ggplot2 0.9.3.1
  require(gridExtra)

  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep('strip_t', names(gg))

  for(ii in seq_along(label.list))  {
    modgrob <- grid::getGrob(gg[[strips[ii]]], 'strip.text', 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- grid::editGrob(modgrob, label=label.list[[ii]])
  }

  g$grobs <- gg
  class(g) <- c('arrange', 'ggplot', class(g)) 
  return(g)
}



pathophys.colors <- c('#005200', '#B20000', '#8F008F', '#0000B2', '#E68A00', 'black')
msig.strip.text <- 'Gene\u2014{MSigDB Collection}\u2014Gene\u2014Disease DWPC'

# AUROC plot helpers
PerfPlot <- function(gg, ymin.auroc) {
  gg <- SetGGTheme(gg) +
  facet_grid(. ~ panel, scales='free_x', space='free_x') +
  geom_hline(aes(yintercept=0.5), color='grey', linetype='solid') +
  stat_summary(fun.y='MeanConfInt', geom='line', color='grey', size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  xlab(NULL) + ylab('AUROC') +
  scale_y_continuous(limits=c(ymin.auroc, 1), breaks=seq(0, 1, .2), expand=c(0.03, 0)) +
  scale_colour_manual(values=pathophys.colors, name='Pathophysiology')
  return(gg)
}

NAtoFALSE <- function(vec) {
  vec[is.na(vec)] <- FALSE
  return(vec)
}

AddPanelColumn <- function(auroc.df) {
  auroc.df$panel <- 'Model'
  auroc.df[NAtoFALSE(auroc.df$metric == 'Path Count'), 'panel'] <- 'Path Count'
  auroc.df[NAtoFALSE(auroc.df$metric == 'DWPC (w=0.4)'), 'panel'] <- 'Degree-Weighted Path Count'
  is.msig.auc <- NAtoFALSE(substr(auroc.df$name, 1, 1) == '{')
  auroc.df[is.msig.auc, 'panel'] <- msig.strip.text
  return(auroc.df)
}


MeanConfInt <- function(x) {t.test(x)$conf.int[1:2]}

PlotAUROCs <- function(auroc.df, perm.df) {
  # Plots to graphics device
  auroc.df <- AddPanelColumn(auroc.df)
  ymin.auroc <- min(subset(auroc.df, auroc != 0)$auroc)
  emdash.size <- 8
  emdash.perm.size <- 6
  jitter.width <- 0.3
  point.size <- 1.75


  if (! missing(perm.df)) {
    perm.df <- AddPanelColumn(perm.df)
    auroc.df$permuted <- FALSE
    perm.df$permuted <- TRUE
		auroc.df <- rbind(auroc.df, subset(perm.df, breadth == 'global'))
  }


  ## MSigDB Plot
  msig.df <- subset(auroc.df, panel == msig.strip.text | panel == 'Model')
  msig.df$name <- gsub('[{}]', '', msig.df$name)
  msig.disease.df <- subset(msig.df, breadth == 'disease_specific')
  msig.global.df <- subset(msig.df, breadth == 'global')
  # Order by global AUROC
  if (! missing(perm.df)) {
		msig.global.perm.df <- subset(msig.global.df, permuted)
		msig.global.df <- subset(msig.global.df, ! permuted)
	}
  msig.levels <- as.character(msig.global.df$name[order(msig.global.df$auroc)])
  msig.disease.df$name <- factor(msig.disease.df$name, levels=msig.levels)
  msig.global.df$name <- factor(msig.global.df$name, levels=msig.levels)
  if (! missing(perm.df)) {
		msig.global.perm.df$name <- factor(msig.global.perm.df$name, levels=msig.levels)
	}

  set.seed(0); msig.plot <- ggplot(msig.disease.df, aes(name, auroc)) 
  msig.plot <- PerfPlot(msig.plot, ymin.auroc) +
    geom_point(data=msig.global.df, size=emdash.size, shape='\u2014', color=Solar('base02')) + 
    geom_point(aes(color=disease_pathophys), position=position_jitter(width=jitter.width), 
      alpha=0.7, size=point.size, show_guide=FALSE)
  if (! missing(perm.df)) {
    msig.plot <- msig.plot + geom_point(data=msig.global.perm.df, 
      size=emdash.perm.size, shape='\u2014', color=Solar('violet'))
  }

  ## NonMSigDB Plot
  nonmsig.df <- subset(auroc.df, panel != msig.strip.text)
  nonmsig.df$panel <- factor(nonmsig.df$panel, levels=c('Degree-Weighted Path Count', 'Path Count', 'Model'))
  nonmsig.disease.df <- subset(nonmsig.df, breadth == 'disease_specific' & panel != 'PCt')
  nonmsig.global.df <- subset(nonmsig.df, breadth == 'global')
  # Order by global AUROC
  if (! missing(perm.df)) {
		nonmsig.global.perm.df <- subset(nonmsig.global.df, permuted)
		nonmsig.global.df <- subset(nonmsig.global.df, ! permuted)
	}
  nonmsig.levels <- unique(as.character(nonmsig.global.df$name[order(nonmsig.global.df$auroc)]))
  nonmsig.disease.df$name <- factor(nonmsig.disease.df$name, levels=nonmsig.levels)
  nonmsig.global.df$name <- factor(nonmsig.global.df$name, levels=nonmsig.levels)
  if (! missing(perm.df)) {
		nonmsig.global.perm.df$name <- factor(nonmsig.global.perm.df$name, levels=nonmsig.levels)
	}


  set.seed(0); nonmsig.plot <- ggplot(nonmsig.disease.df, aes(name, auroc))
  nonmsig.plot <- PerfPlot(nonmsig.plot, ymin.auroc) +
    geom_point(data=nonmsig.global.df, size=emdash.size, shape='\u2014', color=Solar('base02')) + 
    geom_point(aes(color=disease_pathophys), position=position_jitter(width=jitter.width),
      alpha=0.7, size=point.size) + 
    theme(legend.key=element_rect(linetype='blank')) +
    theme(legend.margin=grid::unit(1, 'points'))
  if (! missing(perm.df)) {
    nonmsig.plot <- nonmsig.plot + geom_point(data=nonmsig.global.perm.df, 
      size=emdash.perm.size, shape='\u2014', color=Solar('violet'))
  }

  gg = gridExtra::grid.arrange(msig.plot, nonmsig.plot, ncol=1, widths=c(1, 1))
  return(gg)
}

