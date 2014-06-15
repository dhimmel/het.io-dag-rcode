#http://ije.oxfordjournals.org/content/suppl/2012/02/15/dyr241.DC1/appendix.pdf
# OR to RR http://www.r-bloggers.com/how-to-convert-odds-ratios-to-relative-risks/




BayesianFDP <- function(log.odds, variance, prior.variance, prior) {
  # This function calculates BFDP, the approximate Pr( H0 | thetahat ),
  # given an estiamte of the log relative risk, thetahat, the variance of
  # this estimate, V, the prior variance, W, and the prior probability of
  # a non-null association.
  # http://faculty.washington.edu/jonno/BFDP.R
  # Bayesian false-discovery probability
  # http://faculty.washington.edu/jonno/software.html
  # doi: 10.1093/ije/dyr241
  # http://ije.oxfordjournals.org/content/suppl/2012/02/15/dyr241.DC1/appendix.pdf
  H0.prob <- dnorm(log.odds, m=0, s=sqrt(variance))
  post.var <- variance + prior.variance
  H1.prob <- dnorm(log.odds, m=0, s=sqrt(post.var))
  bayes.factor <- H0.prob / H1.prob
  prior.odds.null <- (1 - prior) / prior
  bfdp <- bayes.factor * prior.odds.null / (bayes.factor * prior.odds.null + 1)
  list('bayes.factor'=bayes.factor, 'H0.prob'=H0.prob, 'H1.prob'=H1.prob, 'bfdp'=bfdp)
}

PriorVarianceFromSampleSize <- function(sample.size, minor.allele.freq) {
  # sample.size=5000 refers to 5000 cases and 5000 controls
  variance <- 1 / (sample.size * minor.allele.freq * (1 - minor.allele.freq))
  return(variance)
}

BayesianAssociation <- function(p.value, odds.ratio, minor.allele.freq, sample.size, snps) {
  log.odds <- log(odds.ratio)
  prior <- 1 / snps
  zscore <- qnorm(1 - p.value / 2) # Calculate the z score
  log.odds.variance <- (log.odds / zscore) ^ 2 # and the asymptotic variance
  prior.variance <- PriorVarianceFromSampleSize(sample.size, minor.allele.freq)
  bayesian.fdp <- BayesianFDP(log.odds, log.odds.variance, prior.variance, prior)
  return(bayesian.fdp)
}







BayesianAssociation(p.value=3e-10, odds.ratio=1.34, minor.allele.freq=0.27, sample.size=45, snps=2.5e-5)




PriorVarianceFromCI <- function() {
  # requires realtive risk
  variance <- ((log(RRhi) - log(RRhat)) / 1.96) ^ 2
  return(variance)
}

















WABF <- function(theta.hat, theta.hat.sd, prior.sd) {
  # Returns the Wakefield Approximate Bayes Factor for a GWAS association.
  # Input the estimated log odds ratio (theta.hat), the standard deviation of that 
  # estimate (theta.hat.sd) and the standard deviation of the prior distribution
  # of effect sizes.
  # From doi:10.1038/nrg2615 (page 683, formula 6)
  Z <- theta.hat / theta.hat.sd
  theta.hat.var <- theta.hat.sd ^ 2
  prior.var <- prior.sd ^ 2
  total.var <- theta.hat.var + prior.var
  wabf <- sqrt(theta.hat.var / (total.var)) * exp(prior.var * Z ^ 2 / (2 * total.var))
  return(wabf)
}


#ThetaHatSDfromP <- function(odds.ratio, p.value) {
#  zscore <- qnorm(1 - p.value / 2)
#  theta.hat.sd <- (log(odds.ratio) / zscore)
#  return(theta.hat.sd)
#}


ThetaHatSDfromP <- function(odds.ratio, mlog10.pval) {
  log10.pval <- -mlog10.pval
  log.pval <- log10.pval / log10(exp(1))
  zscore <- abs(qnorm(log.pval - log(2), log.p=TRUE))
  theta.hat.sd <- (log(odds.ratio) / zscore)
  return(theta.hat.sd)
}



ThetaHatSDfromCI <- function(odds.ratio, bound, level=0.95) {
  log.odds.ratio <- log(odds.ratio)
  log.bound <- log(bound)
  span <- abs(log.bound - log.odds.ratio)
  theta.hat.sd <- span / 1.96
  return(theta.hat.sd)
}

CalculatePPA <- function(bayes.factor, prior.prob) {
  posterior.odds <- bayes.factor * prior.prob / (1 - prior.prob)
  ppa <- posterior.odds / (1 + posterior.odds)
  return(ppa)
}

PosteriorProbability <- function(odds.ratio, mlog10.pval, theta.prior.sd=0.2, prior.prob=1e-6) {
  # posterior probability of association
  theta.hat <- log(odds.ratio)
  theta.hat.sd <- ThetaHatSDfromP(odds.ratio, mlog10.pval)
  wabf <- WABF(theta.hat, theta.hat.sd, theta.prior.sd)
  ppa <- CalculatePPA(wabf, prior.prob)
  return(ppa)
}

PosteriorProbability(odds.ratio=1.34, mlog10.pval=10.1, theta.prior.sd=0.2, prior.prob=1e-6)
PosteriorProbability(odds.ratio=1.64, p.value=2e-18, theta.prior.sd=0.2, prior.prob=1e-6)


PosteriorProbability(odds.ratio=1.09, mlog10.pval=-log10(1e-6), theta.prior.sd=0.2, prior.prob=1e-5)


gcat.df <- read.delim('/home/dhimmels/Documents/serg/data-sources/gwas-catalog/140205/bayesian/catalog-rows.txt', stringsAsFactors=FALSE, check.names=FALSE)

library(ggplot2)

ggplot(gcat.df, aes(mlog_pval, ppa, size=cases)) +
  geom_point() + xlim(c(5, 15)) +
  theme_bw()

ggplot(gcat.df, aes(mlog_pval, ppa, color=merged_sample_size > 1000)) +
  geom_point() + xlim(c(5, 15)) +
  theme_bw()






