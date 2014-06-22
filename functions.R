
# Increase terminal width in interactive sessions
terminal.width <- Sys.getenv('COLUMNS')
terminal.width <- ifelse(terminal.width == '', 80, terminal.width)
options(width=terminal.width)

# plotting parameters
width.full <- 6.83
width.half <- 3.27

breaks.roc <- seq(0, 1, by=0.2)
breaks.prc <- seq(0, 1, by=0.2)

# Colors
strip.fill <- '#fdf6e3' # solarized website background


'#300A24' # Blackberry (ubuntu terminal)
'#F8F6F7' # lighter blackberry


ReadFeatures <- function(feature.dir) {
  # Read and combine the feature files in feature.dir with names beginning with 'DOID_'.
  feature.filenames <- list.files(feature.dir, pattern='DOID_')
  feature.df.list <- list()
  for (feature.filename in feature.filenames) {
    doid_code <- gsub('.txt.gz', '', feature.filename)
    cat(sprintf('Reading features for %s\n', doid_code))
    feature.path <- file.path(feature.dir, feature.filename)
    feature.df <- read.delim(feature.path, check.names=FALSE)
    feature.df.list[[doid_code]] <- feature.df
  }
  return(do.call(rbind, feature.df.list))
}


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

GLMNetCoef <- function(cv.glm, X, y, prepend='', ...) {
  lambda <- cv.glm$lambda.1se
  coef.vec <- coef(cv.glm, s=lambda)[, 1]
  zcoef.vec <- c(0, coef.vec[-1] * apply(X, 2, sd) / sd(y))
  coef.df <- data.frame('feature'=c('intercept', colnames(X)), ...,
    'coef'=coef.vec, 'zcoef'=zcoef.vec, row.names=NULL)
  colnum <- ncol(coef.df)
  colnames(coef.df)[c(colnum - 1, colnum)] <- paste0(prepend, c('coef', 'zcoef'))
  return(coef.df)
}

## Plotting functions
font.dir <- '/home/dhimmels/Documents/serg/gene-disease-hetnet/fonts'

arial <- c(file.path(font.dir, 'Arial.afm'  ),
           file.path(font.dir, 'Arial_Bold.afm'),
           file.path(font.dir, 'Arial_Italic.afm' ),
           file.path(font.dir, 'Arial_Bold_Italic.afm'))

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
  scale_y_continuous(breaks=breaks.prc, expand=c(0.03, 0)) +
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

