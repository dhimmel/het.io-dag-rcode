



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

## Plotting functions

ChrRound <- function(x, digits=2) {
  # Round x to digits and return as a character with trailing zeros.
  x.rounded <- round(x, digits)
  sprintf.chr <- sprintf('%%.%sf', digits)
  x.chr <- sprintf(sprintf.chr, x.rounded)
  return(x.chr)
}


SetGGTheme <- function(gg) {
  gg <- gg + theme_bw() + theme(legend.title=element_text(face='bold.italic'))
  gg <- gg + theme(plot.margin=grid::unit(c(2, 2, 2, 2), 'points'))
  return(gg)
}


ggROC <- function(gg) {
  gg <- SetGGTheme(gg) + 
    geom_abline(slope=1, color='grey') +
    theme(legend.background=element_rect(color='grey60', size=0.25)) +
    theme(legend.key.width=grid::unit(1.5, 'lines')) +
    theme(legend.key.height=grid::unit(0.95, 'lines')) +
    theme(legend.key=element_rect(linetype='blank')) +
    theme(legend.justification=c(1,0), legend.position=c(1,0)) +
    xlab('False Positive Rate') + ylab('Recall') + coord_fixed()
  return(gg)
}

