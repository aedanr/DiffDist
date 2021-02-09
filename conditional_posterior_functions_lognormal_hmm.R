# Per-gene conditional log posterior mean function
# Takes vectors of per-gene sample means and mean and dispersion estimates
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


# Per-gene conditional log posterior dispersion function
# Takes vector of counts for each gene to calculate sum over log gamma for each gene, 
# and vectors of per-gene sample means and mean and dispersion estimates
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
  lgammasum <- colSums(lgamma(counts + rep(x, each = nrow(counts))))
  out[index] <- 
    lgammasum[index] - 
    n * lgamma(x[index]) - 
    n * (sample.means[index] + x[index]) * log(1 + means[index]*disps[index]) + 
    (n*sample.means[index] - 1) * log(disps[index]) - 
    (log(disps[index]) - prior.location)^2 / (2*prior.scale)
  return(out)
}


# Function for sd of normal conditional posteriors for lognormal location parameters
# Takes vector of posterior indicators
prior.location.posterior.sd <- function(hyperprior.var, 
                                        prior.scale, 
                                        z) {
  return(sqrt(1 / (1/hyperprior.var + 
                     sum(1 + z) / prior.scale)))
}


# Function for mean of normal conditional posteriors for lognormal location parameters
# Takes vectors of mean or dispersion estimates and posterior indicators
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


# Conditional log posterior function for scale parameter of prior on mean
# Takes vectors of per-gene mean estimates and mixture components
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


# Conditional log posterior function for scale parameter of prior on dispersion
# Takes vectors of per-gene dispersion estimates and mixture components
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


# Function for log posterior probability z=0
# Takes vector of counts for each gene to calculate sum over log gamma for each gene, 
# and vectors of per-gene sample means and mean and dispersion estimates
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
  lgammasum <- colSums(lgamma(counts + rep(x, each = nrow(counts))))
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


# Function for log posterior probability z=1
# Takes vector of counts for each gene to calculate sums over log gamma for each gene 
# for each group, and vectors of per-gene sample means and mean and dispersion estimates for each group
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
  lgammasum1 <- colSums(lgamma(counts1 + rep(x1, each = nrow(counts1))))
  lgammasum2 <- colSums(lgamma(counts2 + rep(x2, each = nrow(counts2))))
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


# Function to calculate probabilities for posterior Bernoulli distributions for 
# mixture components (exponentiates and normalises calculated probabilities)
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
