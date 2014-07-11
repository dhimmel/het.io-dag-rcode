
CROC <- function(score, status, directory, alpha=7) {
  # Requires CROC and sympy python package installations
  # Online documentation: http://swami.wustl.edu/CROC/
  # Python package https://pypi.python.org/pypi/CROC/
  # Source publication doi: 10.1093/bioinformatics/btq140
  
  # Establish temporary directory
  dir0 <- getwd()
  if (missing(directory)) {
    directory <- file.path(dir0, Sys.getpid())}

  created.dir <- dir.create(directory)
  if (! created.dir) {
    stop(paste0('preexisting CROC temp dir: '), directory)}
  setwd(directory)

  # Establish file paths
  files <- list(directory=directory, scores='input.scored-label',
    output='output.curve', random='random.curve', info='info.txt')

  # Write input
  scored.label.df <- data.frame('score'=score, 'label'=status)
  write.table(scored.label.df, files$scores, row.names=FALSE, col.names=FALSE, sep='\t')
  
  # Create system call to execute CROC python analysis
  sys2.args <- c(
    '--tie_mode', 'smooth',
    '--curve_type', 'roc',
    '--transform', shQuote(sprintf('Exponential(%s)', alpha))
  )
  sys2.out <- system2('croc-curve', args=sys2.args, stdin=files$scores,
    stdout=files$output, stderr=files$info)

  # Read output
  info <- scan(files$info, what='character', sep='\n')
  auc <- as.numeric(sub('Area Under Curve =[ ]+', '', info[1]))
  curve.df <- read.table(files$output, col.names=c('x', 'y'))
  random.df <- data.frame('x'=seq(0, 1, 0.005))
  if (alpha == 0) {
    random.df$y <- random.df$x
  } else {
    # formula 5 from paper
    random.df$y <- -log(1 - random.df$x * (1 - exp(-alpha))) / alpha
  }
  
  # Remove temp dir and reset working directory
  unlink(directory, recursive=TRUE)
  setwd(dir0)

  # Create list to return
  croc <- list()
  croc$args <- sys2.args
  croc$auc <- auc
  croc$curve.df <- curve.df
  croc$random.df <- random.df
  return(croc)
}

ExpTransform <- function(fpr, alpha) {
  (1 - exp(-alpha * fpr)) / (1 - exp(-alpha))
}

FindAlpha <- function(fpr0, fpr1) {
  minimize <- function(alpha, fpr0, fpr1) {(ExpTransform(fpr0, alpha) - fpr1) ^ 2}
  alpha <- optimize(minimize, fpr0=fpr0, fpr1=fpr1, lower=0, upper=750, tol=1e-20)$minimum
  return(alpha)
}

