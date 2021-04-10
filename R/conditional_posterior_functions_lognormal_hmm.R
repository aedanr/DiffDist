#' Per-gene conditional log posterior mean function.
#'
#' Called by MCMC functions.
#' @param n Number of samples
#' @param sample.means Vector of per-gene sample means
#' @param means Vector of current per-gene mean estimates
#' @param disps Vector of current per-gene dispersion estimates
#' @param prior.location Current estimate of prior location parameter
#' @param prior.scale Current estimate of prior scale parameter
#' @return Vector of conditional posterior mean densities on log scale
#' @export
mean.conditional.log.posterior <- function(n,
                                           sample.means,
                                           means,
                                           disps,
                                           prior.location,
                                           prior.scale) {
  x <- 1/disps
  out <- rep(-1e20, length(means))
  index <- which(means>0)
  out[index] <-
    -n *
    (sample.means[index] + x[index]) *
    log(1 + means[index]*disps[index]) +
    (n*sample.means[index] - 1) *
    log(means[index]) -
    (log(means[index]) - prior.location)^2 / (2*prior.scale)
  return(out)
}


#' Per-gene conditional log posterior dispersion function. Called by MCMC functions.
#' @param genes Number of genes
#' @param counts Matrix of counts with genes in rows and samples in columns
#' @param n Number of samples
#' @param sample.means Vector of per-gene sample means
#' @param means Vector of current per-gene mean estimates
#' @param disps Vector of current per-gene dispersion estimates
#' @param prior.location Current estimate of prior location parameter
#' @param prior.scale Current estimate of prior scale parameter
#' @return Vector of conditional posterior dispersion densities on log scale
#' @export
disp.conditional.log.posterior = function(genes,
                                          counts,
                                          n,
                                          sample.means,
                                          means,
                                          disps,
                                          prior.location,
                                          prior.scale) {
  x <- 1/disps
  out <- rep(-1e20,genes)
  index <- which(disps>0)
  lgammasum <- rowSums(lgamma(counts + rep(x, ncol(counts))))
  out[index] <-
    lgammasum[index] -
    n * lgamma(x[index]) -
    n * (sample.means[index] + x[index]) * log(1 + means[index]*disps[index]) +
    (n*sample.means[index] - 1) * log(disps[index]) -
    (log(disps[index]) - prior.location)^2 / (2*prior.scale)
  return(out)
}


#' Function to compute standard deviation of conditional posteriors for prior
#' location parameters. Used for Gibbs sampling for prior location parameters.
#' @param hyperprior.var Variance of hyperprior on prior location parameter
#' @param prior.scale Current estimate of prior scale parameter
#' @param z Vector of current per-gene mixture component indicator variables
#' @return Standard deviation of conditional posterior distribution of the given prior
#' location parameter
#' @export
prior.location.posterior.sd <- function(hyperprior.var,
                                        prior.scale,
                                        z) {
  return(sqrt(1 / (1/hyperprior.var +
                     sum(1 + z) / prior.scale)))
}


#' Function to compute mean of conditional posteriors for prior location parameters.
#' Used for Gibbs sampling for prior location parameters.
#' @param hyperprior.mean Mean of hyperprior on prior location parameter
#' @param hyperprior.var Variance of hyperprior on prior location parameter
#' @param prior.scale Current estimate of prior scale parameter
#' @param parameters0 Vector of current overall mean or dispersion estimates
#' @param parameters1 Vector of current group 1 mean or dispersion estimates
#' @param parameters2 Vector of current group 2 mean or dispersion estimates
#' @param z Vector of current per-gene mixture component indicator variables
#' @return Mean of conditional posterior distribution of the given prior location
#' parameter
#' @export
prior.location.posterior.mean <- function(hyperprior.mean,
                                          hyperprior.var,
                                          prior.scale,
                                          parameters0,
                                          parameters1,
                                          parameters2,
                                          z) {
  return(
    (hyperprior.mean / hyperprior.var +
       sum((1-z)*log(parameters0) + z*(log(parameters1) + log(parameters2))) /
       prior.scale) /
      (1/hyperprior.var + sum(1+z)/prior.scale)
  )
}


#' Conditional log posterior function for scale parameter of mean prior. Called by
#' MCMC functions.
#' @param prior.scale Current estimate of prior scale parameter
#' @param prior.location Current estimate of prior location parameter
#' @param g Number of genes
#' @param means0 Vector of current overall per-gene mean estimates
#' @param means1 Vector of current group 1 per-gene mean estimates
#' @param means2 Vector of current group 2 per-gene mean estimates
#' @param z Vector of current per-gene mixture component indicator variables
#' @return Conditional posterior density of scale parameter of mean prior, on log scale
#' @export
mean.prior.scale.log.posterior <- function(prior.scale,
                                           prior.location,
                                           g,
                                           means0,
                                           means1,
                                           means2,
                                           z) {
  if (prior.scale <= 0) {
    return(-1e20)
  } else {
    return(
      (2 - (g + sum(z))/2) * log(prior.scale) -
        0.8 * prior.scale - sum(
          (1-z)*(log(means0) - prior.location)^2 +
            z*((log(means1) - prior.location)^2 +
                 (log(means2) - prior.location)^2)
        ) / (2*prior.scale)
    )
  }
}


#' Conditional log posterior function for scale parameter of dispersion prior. Called
#' by MCMC functions.
#' @param prior.scale Current estimate of prior scale parameter
#' @param prior.location Current estimate of prior location parameter
#' @param g Number of genes
#' @param disps0 Vector of current overall per-gene dispersion estimates
#' @param disps1 Vector of current group 1 per-gene dispersion estimates
#' @param disps2 Vector of current group 2 per-gene dispersion estimates
#' @param z Vector of current per-gene mixture component indicator variables
#' @return Conditional posterior density of scale parameter of mean dispersion, on log
#' scale
#' @export
disp.prior.scale.log.posterior <- function(prior.scale,
                                           prior.location,
                                           g,
                                           disps0,
                                           disps1,
                                           disps2,
                                           z) {
  if (prior.scale <= 0) {
    return(-1e20)
  } else {
    return(
      (4 - (g + sum(z))/2) * log(prior.scale) -
        2 * prior.scale - sum(
          (1-z)*(log(disps0) - prior.location)^2 +
            z*((log(disps1) - prior.location)^2 +
                 (log(disps2) - prior.location)^2)
        ) / (2*prior.scale)
    )
  }
}


#' Function for raw conditional log posterior probability of no differential
#' distribution. Called by \code{posterior.indicator.probabilities()}.
#' @param genes Number of genes
#' @param counts Matrix of counts with genes in rows and samples in columns
#' @param n Number of samples
#' @param sample.means Vector of per-gene sample means
#' @param means Vector of current per-gene mean estimates
#' @param disps Vector of current per-gene dispersion estimates
#' @param mean.prior.location Current estimate of mean prior location parameter
#' @param mean.prior.scale Current estimate of mean prior scale parameter
#' @param disp.prior.location Current estimate of dispersion prior location parameter
#' @param disp.prior.scale Current estimate of dispersion prior scale parameter
#' @param lambda Current estimate of proportion of differentially distributed genes
#' @return Vector of un-normalised per-gene conditional posterior probabilities of no
#' differential distribution, on log scale
#' @export
pz0 <- function(genes,
                counts,
                n,
                sample.means,
                means,
                disps,
                mean.prior.location,
                mean.prior.scale,
                disp.prior.location,
                disp.prior.scale,
                lambda) {
  x <- 1/disps
  lgammasum <- rowSums(lgamma(counts + rep(x, ncol(counts))))
  return(lgammasum -
           n * lgamma(x) -
           n * (sample.means + x) * log(1 + means * disps) +
           (n * sample.means - 1) * log(means * disps) +
           log(2*pi) +
           log(mean.prior.scale)/2 +
           log(disp.prior.scale)/2 -
           (log(means) - mean.prior.location)^2 / (2*mean.prior.scale) -
           (log(disps) - disp.prior.location)^2 / (2*disp.prior.scale) +
           log(1-lambda))
}


#' Function for raw conditional log posterior probability of differential
#' distribution. Called by \code{posterior.indicator.probabilities()}.
#' @param genes Number of genes
#' @param counts1 Matrix of counts in group 1 with genes in rows and samples in columns
#' @param counts2 Matrix of counts in group 2
#' @param n1 Number of samples in group 1
#' @param n2 Number of samples in group 2
#' @param sample.means1 Vector of per-gene group 1 sample means
#' @param sample.means2 Vector of per-gene group 2 sample means
#' @param means1 Vector of current per-gene group 1 mean estimates
#' @param means2 Vector of current per-gene group 2 mean estimates
#' @param disps1 Vector of current per-gene group 1 dispersion estimates
#' @param disps2 Vector of current per-gene group 2 dispersion estimates
#' @param mean.prior.location Current estimate of mean prior location parameter
#' @param mean.prior.scale Current estimate of mean prior scale parameter
#' @param disp.prior.location Current estimate of dispersion prior location parameter
#' @param disp.prior.scale Current estimate of dispersion prior scale parameter
#' @param lambda Current estimate of proportion of differentially distributed genes
#' @return Vector of un-normalised per-gene conditional posterior probabilities of
#' differential distribution, on log scale
#' @export
pz1 = function(genes,
               counts1,
               counts2,
               n1,
               n2,
               sample.means1,
               sample.means2,
               means1,
               means2,
               disps1,
               disps2,
               mean.prior.location,
               mean.prior.scale,
               disp.prior.location,
               disp.prior.scale,
               lambda) {
  x1 <- 1/disps1
  x2 <- 1/disps2
  lgammasum1 <- rowSums(lgamma(counts1 + rep(x1, ncol(counts1))))
  lgammasum2 <- rowSums(lgamma(counts2 + rep(x2, ncol(counts2))))
  return(lgammasum1 -
           n1 * lgamma(x1) -
           n1 * (sample.means1 + x1) * log(1 + means1 * disps1) +
           (n1 * sample.means1 - 1) * log(means1 * disps1) +
           lgammasum2 -
           n2 * lgamma(x2) -
           n2 * (sample.means2 + x2) * log(1 + means2 * disps2) +
           (n2 * sample.means2 - 1) * log(means2 * disps2) -
           ((log(means1) - mean.prior.location)^2 +
              (log(means2) - mean.prior.location)^2) / (2*mean.prior.scale) -
           ((log(disps1) - disp.prior.location)^2 +
              (log(disps2) - disp.prior.location)^2) / (2*disp.prior.scale) +
           log(lambda))
}


#' Function for posterior probability of differential distribution. Called by MCMC
#' functions.
#' @param genes Number of genes
#' @param counts0 Overall count matrix with genes in rows and samples in columns
#' @param counts1 Matrix of counts in group 1
#' @param counts2 Matrix of counts in group 2
#' @param n Total number of samples
#' @param n1 Number of samples in group 1
#' @param n2 Number of samples in group 2
#' @param sample.means0 Vector of overall sample means
#' @param sample.means1 Vector of per-gene group 1 sample means
#' @param sample.means2 Vector of per-gene group 2 sample means
#' @param means0 Vector of current overall mean estimates
#' @param means1 Vector of current per-gene group 1 mean estimates
#' @param means2 Vector of current per-gene group 2 mean estimates
#' @param disps0 Vector of current overall dispersion estimates
#' @param disps1 Vector of current per-gene group 1 dispersion estimates
#' @param disps2 Vector of current per-gene group 2 dispersion estimates
#' @param mean.prior.location Current estimate of mean prior location parameter
#' @param mean.prior.scale Current estimate of mean prior scale parameter
#' @param disp.prior.location Current estimate of dispersion prior location parameter
#' @param disp.prior.scale Current estimate of dispersion prior scale parameter
#' @param lambda Current estimate of proportion of differentially distributed genes
#' @return Vector of posterior probabilities of differential distribution
#' @export
posterior.indicator.probabilities <- function(genes,
                                              counts0,
                                              counts1,
                                              counts2,
                                              n,
                                              n1,
                                              n2,
                                              sample.means0,
                                              sample.means1,
                                              sample.means2,
                                              means0,
                                              means1,
                                              means2,
                                              disps0,
                                              disps1,
                                              disps2,
                                              mean.prior.location,
                                              mean.prior.scale,
                                              disp.prior.location,
                                              disp.prior.scale,
                                              lambda) {
  return(1 / (1 + exp(pz0(genes=genes,
                          counts=counts0,
                          n=n,
                          sample.means=sample.means0,
                          means=means0,
                          disps=disps0,
                          mean.prior.location=mean.prior.location,
                          mean.prior.scale=mean.prior.scale,
                          disp.prior.location=disp.prior.location,
                          disp.prior.scale=disp.prior.scale,
                          lambda=lambda) -
                        pz1(genes=genes,
                            counts1=counts1,
                            counts2=counts2,
                            n1=n1,
                            n2=n2,
                            sample.means1=sample.means1,
                            sample.means2=sample.means2,
                            means1=means1,
                            means2=means2,
                            disps1=disps1,
                            disps2=disps2,
                            mean.prior.location=mean.prior.location,
                            mean.prior.scale=mean.prior.scale,
                            disp.prior.location=disp.prior.location,
                            disp.prior.scale=disp.prior.scale,
                            lambda=lambda)
  )
  )
  )
}
