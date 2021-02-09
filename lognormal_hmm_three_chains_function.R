ln_hmm_3_chains <- function(counts, 
                            groups, 
                            chain.length, 
                            thin=1, 
                            inits1, 
                            inits2, 
                            inits3, 
                            mean.proposal.scales0=rep(0.2, ncol(counts)), 
                            mean.proposal.scales1=rep(0.2, ncol(counts)), 
                            mean.proposal.scales2=rep(0.2, ncol(counts)), 
                            disp.proposal.scales0=rep(0.5, ncol(counts)), 
                            disp.proposal.scales1=rep(0.5, ncol(counts)), 
                            disp.proposal.scales2=rep(0.5, ncol(counts)), 
                            mean.prior.scale.proposal.sd=0.1, 
                            disp.prior.scale.proposal.sd=0.4) {
  
  genes <- ncol(counts)
  counts1 <- counts[groups==1,]
  counts2 <- counts[groups==2,]
  samples0 <- nrow(counts)
  samples1 <- nrow(counts1)
  samples2 <- nrow(counts2)
  sample.means0 <- colMeans(counts)
  sample.means1 <- colMeans(counts1)
  sample.means2 <- colMeans(counts2)
  sqrt.mean.proposal.scales0 <- sqrt(mean.proposal.scales0)
  sqrt.mean.proposal.scales1 <- sqrt(mean.proposal.scales1)
  sqrt.mean.proposal.scales2 <- sqrt(mean.proposal.scales2)
  sqrt.disp.proposal.scales0 <- sqrt(disp.proposal.scales0)
  sqrt.disp.proposal.scales1 <- sqrt(disp.proposal.scales1)
  sqrt.disp.proposal.scales2 <- sqrt(disp.proposal.scales2)
  
  # Create empty posterior sample matrices and vectors ####
  posterior.means0.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means1.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means2.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means0.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means1.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means2.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means0.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means1.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.means2.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps0.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps1.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps2.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps0.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps1.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps2.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps0.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps1.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.disps2.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.mean.prior.location.1 <- numeric(chain.length/thin)
  posterior.mean.prior.location.2 <- numeric(chain.length/thin)
  posterior.mean.prior.location.3 <- numeric(chain.length/thin)
  posterior.mean.prior.scale.1 <- numeric(chain.length/thin)
  posterior.mean.prior.scale.2 <- numeric(chain.length/thin)
  posterior.mean.prior.scale.3 <- numeric(chain.length/thin)
  posterior.disp.prior.location.1 <- numeric(chain.length/thin)
  posterior.disp.prior.location.2 <- numeric(chain.length/thin)
  posterior.disp.prior.location.3 <- numeric(chain.length/thin)
  posterior.disp.prior.scale.1 <- numeric(chain.length/thin)
  posterior.disp.prior.scale.2 <- numeric(chain.length/thin)
  posterior.disp.prior.scale.3 <- numeric(chain.length/thin)
  posterior.indicators.1 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.indicators.2 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.indicators.3 <- matrix(nrow=chain.length/thin, ncol=genes)
  posterior.proportion.1 <- numeric(chain.length/thin)
  posterior.proportion.2 <- numeric(chain.length/thin)
  posterior.proportion.3 <- numeric(chain.length/thin)
  
  # Initial values ####
  current.posterior.means0.1 <- inits1$means0
  current.posterior.means1.1 <- inits1$means1
  current.posterior.means2.1 <- inits1$means2
  current.posterior.means0.2 <- inits2$means0
  current.posterior.means1.2 <- inits2$means1
  current.posterior.means2.2 <- inits2$means2
  current.posterior.means0.3 <- inits3$means0
  current.posterior.means1.3 <- inits3$means1
  current.posterior.means2.3 <- inits3$means2
  current.posterior.disps0.1 <- inits1$disps0
  current.posterior.disps1.1 <- inits1$disps1
  current.posterior.disps2.1 <- inits1$disps2
  current.posterior.disps0.2 <- inits2$disps0
  current.posterior.disps1.2 <- inits2$disps1
  current.posterior.disps2.2 <- inits2$disps2
  current.posterior.disps0.3 <- inits3$disps0
  current.posterior.disps1.3 <- inits3$disps1
  current.posterior.disps2.3 <- inits3$disps2
  current.posterior.mean.prior.location.1 <- inits1$mean.prior.location
  current.posterior.mean.prior.location.2 <- inits2$mean.prior.location
  current.posterior.mean.prior.location.3 <- inits3$mean.prior.location
  current.posterior.mean.prior.scale.1 <- inits1$mean.prior.scale
  current.posterior.mean.prior.scale.2 <- inits2$mean.prior.scale
  current.posterior.mean.prior.scale.3 <- inits3$mean.prior.scale
  current.posterior.disp.prior.location.1 <- inits1$disp.prior.location
  current.posterior.disp.prior.location.2 <- inits2$disp.prior.location
  current.posterior.disp.prior.location.3 <- inits3$disp.prior.location
  current.posterior.disp.prior.scale.1 <- inits1$disp.prior.scale
  current.posterior.disp.prior.scale.2 <- inits2$disp.prior.scale
  current.posterior.disp.prior.scale.3 <- inits3$disp.prior.scale
  current.posterior.indicators.1 <- numeric(genes)
  current.posterior.indicators.2 <- numeric(genes)
  current.posterior.indicators.3 <- numeric(genes)
  current.posterior.proportion.1 = 0.5
  current.posterior.proportion.2 = 0.5
  current.posterior.proportion.3 = 0.5
  
  # Create acceptance rate vectors/variables ####
  accept.means0.1 <- numeric(genes)
  accept.means1.1 <- numeric(genes)
  accept.means2.1 <- numeric(genes)
  accept.means0.2 <- numeric(genes)
  accept.means1.2 <- numeric(genes)
  accept.means2.2 <- numeric(genes)
  accept.means0.3 <- numeric(genes)
  accept.means1.3 <- numeric(genes)
  accept.means2.3 <- numeric(genes)
  accept.disps0.1 <- numeric(genes)
  accept.disps1.1 <- numeric(genes)
  accept.disps2.1 <- numeric(genes)
  accept.disps0.2 <- numeric(genes)
  accept.disps1.2 <- numeric(genes)
  accept.disps2.2 <- numeric(genes)
  accept.disps0.3 <- numeric(genes)
  accept.disps1.3 <- numeric(genes)
  accept.disps2.3 <- numeric(genes)
  accept.mean.prior.scale.1 <- 0
  accept.mean.prior.scale.2 <- 0
  accept.mean.prior.scale.3 <- 0
  accept.disp.prior.scale.1 <- 0
  accept.disp.prior.scale.2 <- 0
  accept.disp.prior.scale.3 <- 0
  
  # Run MCMC ####
  for (iter in 1:chain.length) {
    # Metropolis updates for overall per-gene means ####
    # Chain 1
    proposed.posterior.means0.1 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.means0.1), 
      sdlog=sqrt.mean.proposal.scales0
    )
    
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(
        n=samples0, 
        sample.means=sample.means0, 
        means=proposed.posterior.means0.1, 
        disps=current.posterior.disps0.1, 
        prior.location=current.posterior.mean.prior.location.1, 
        prior.scale=current.posterior.mean.prior.scale.1
      ) + 
      (log(proposed.posterior.means0.1) + 
         (log(proposed.posterior.means0.1) - log(current.posterior.means0.1))^2 / 
         (2*mean.proposal.scales0)) - 
      mean.conditional.log.posterior(
        n=samples0, 
        sample.means=sample.means0, 
        means=current.posterior.means0.1, 
        disps=current.posterior.disps0.1, 
        prior.location=current.posterior.mean.prior.location.1, 
        prior.scale=current.posterior.mean.prior.scale.1
      ) - 
      (log(current.posterior.means0.1) + 
         (log(current.posterior.means0.1) - log(proposed.posterior.means0.1))^2 / 
         (2*mean.proposal.scales0))
    
    current.posterior.means0.1[replace] <- proposed.posterior.means0.1[replace]
    accept.means0.1[replace] <- accept.means0.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.means0.2 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.means0.2), 
      sdlog=sqrt.mean.proposal.scales0
    )
    
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(
        n=samples0, 
        sample.means=sample.means0, 
        means=proposed.posterior.means0.2, 
        disps=current.posterior.disps0.2, 
        prior.location=current.posterior.mean.prior.location.2, 
        prior.scale=current.posterior.mean.prior.scale.2
      ) + 
      (log(proposed.posterior.means0.2) + 
         (log(proposed.posterior.means0.2) - log(current.posterior.means0.2))^2 / 
         (2*mean.proposal.scales0)) - 
      mean.conditional.log.posterior(
        n=samples0, 
        sample.means=sample.means0, 
        means=current.posterior.means0.2, 
        disps=current.posterior.disps0.2, 
        prior.location=current.posterior.mean.prior.location.2, 
        prior.scale=current.posterior.mean.prior.scale.2
      ) - 
      (log(current.posterior.means0.2) + 
         (log(current.posterior.means0.2) - log(proposed.posterior.means0.2))^2 / 
         (2*mean.proposal.scales0))
    
    current.posterior.means0.2[replace] <- proposed.posterior.means0.2[replace]
    accept.means0.2[replace] <- accept.means0.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.means0.3 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.means0.3), 
      sdlog=sqrt.mean.proposal.scales0
    )
    
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(
        n=samples0, 
        sample.means=sample.means0, 
        means=proposed.posterior.means0.3, 
        disps=current.posterior.disps0.3, 
        prior.location=current.posterior.mean.prior.location.3, 
        prior.scale=current.posterior.mean.prior.scale.3
      ) + 
      (log(proposed.posterior.means0.3) + 
         (log(proposed.posterior.means0.3) - log(current.posterior.means0.3))^2 / 
         (2*mean.proposal.scales0)) - 
      mean.conditional.log.posterior(
        n=samples0, 
        sample.means=sample.means0, 
        means=current.posterior.means0.3, 
        disps=current.posterior.disps0.3, 
        prior.location=current.posterior.mean.prior.location.3, 
        prior.scale=current.posterior.mean.prior.scale.3
      ) - 
      (log(current.posterior.means0.3) + 
         (log(current.posterior.means0.3) - log(proposed.posterior.means0.3))^2 / 
         (2*mean.proposal.scales0))
    
    current.posterior.means0.3[replace] <- proposed.posterior.means0.3[replace]
    accept.means0.3[replace] <- accept.means0.3[replace] + 1/chain.length
    
    # Metropolis updates for group 1 per-gene means ####
    # Chain 1
    proposed.posterior.means1.1 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.means1.1), 
      sdlog=sqrt.mean.proposal.scales1
    )
    
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(
        n=samples1, 
        sample.means=sample.means1, 
        means=proposed.posterior.means1.1, 
        disps=current.posterior.disps1.1, 
        prior.location=current.posterior.mean.prior.location.1, 
        prior.scale=current.posterior.mean.prior.scale.1
      ) + 
      (log(proposed.posterior.means1.1) + 
         (log(proposed.posterior.means1.1) - log(current.posterior.means1.1))^2 / 
         (2*mean.proposal.scales1)) - 
      mean.conditional.log.posterior(
        n=samples1, 
        sample.means=sample.means1, 
        means=current.posterior.means1.1, 
        disps=current.posterior.disps1.1, 
        prior.location=current.posterior.mean.prior.location.1, 
        prior.scale=current.posterior.mean.prior.scale.1
      ) - 
      (log(current.posterior.means1.1) + 
         (log(current.posterior.means1.1) - log(proposed.posterior.means1.1))^2 / 
         (2*mean.proposal.scales1))
    
    current.posterior.means1.1[replace] <- proposed.posterior.means1.1[replace]
    accept.means1.1[replace] <- accept.means1.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.means1.2 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.means1.2), 
      sdlog=sqrt.mean.proposal.scales1
    )
    
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(
        n=samples1, 
        sample.means=sample.means1, 
        means=proposed.posterior.means1.2, 
        disps=current.posterior.disps1.2, 
        prior.location=current.posterior.mean.prior.location.2, 
        prior.scale=current.posterior.mean.prior.scale.2
      ) + 
      (log(proposed.posterior.means1.2) + 
         (log(proposed.posterior.means1.2) - log(current.posterior.means1.2))^2 / 
         (2*mean.proposal.scales1)) - 
      mean.conditional.log.posterior(
        n=samples1, 
        sample.means=sample.means1, 
        means=current.posterior.means1.2, 
        disps=current.posterior.disps1.2, 
        prior.location=current.posterior.mean.prior.location.2, 
        prior.scale=current.posterior.mean.prior.scale.2
      ) - 
      (log(current.posterior.means1.2) + 
         (log(current.posterior.means1.2) - log(proposed.posterior.means1.2))^2 / 
         (2*mean.proposal.scales1))
    
    current.posterior.means1.2[replace] <- proposed.posterior.means1.2[replace]
    accept.means1.2[replace] <- accept.means1.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.means1.3 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.means1.3), 
      sdlog=sqrt.mean.proposal.scales1
    )
    
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(
        n=samples1, 
        sample.means=sample.means1, 
        means=proposed.posterior.means1.3, 
        disps=current.posterior.disps1.3, 
        prior.location=current.posterior.mean.prior.location.3, 
        prior.scale=current.posterior.mean.prior.scale.3
      ) + 
      (log(proposed.posterior.means1.3) + 
         (log(proposed.posterior.means1.3) - log(current.posterior.means1.3))^2 / 
         (2*mean.proposal.scales1)) - 
      mean.conditional.log.posterior(
        n=samples1, 
        sample.means=sample.means1, 
        means=current.posterior.means1.3, 
        disps=current.posterior.disps1.3, 
        prior.location=current.posterior.mean.prior.location.3, 
        prior.scale=current.posterior.mean.prior.scale.3
      ) - 
      (log(current.posterior.means1.3) + 
         (log(current.posterior.means1.3) - log(proposed.posterior.means1.3))^2 / 
         (2*mean.proposal.scales1))
    
    current.posterior.means1.3[replace] <- proposed.posterior.means1.3[replace]
    accept.means1.3[replace] <- accept.means1.3[replace] + 1/chain.length
    
    # Metropolis updates for group 2 per-gene means ####
    # Chain 1
    proposed.posterior.means2.1 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.means2.1), 
      sdlog=sqrt.mean.proposal.scales2
    )
    
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(
        n=samples2, 
        sample.means=sample.means2, 
        means=proposed.posterior.means2.1, 
        disps=current.posterior.disps2.1, 
        prior.location=current.posterior.mean.prior.location.1, 
        prior.scale=current.posterior.mean.prior.scale.1
      ) + 
      (log(proposed.posterior.means2.1) + 
         (log(proposed.posterior.means2.1) - log(current.posterior.means2.1))^2 / 
         (2*mean.proposal.scales2)) - 
      mean.conditional.log.posterior(
        n=samples2, 
        sample.means=sample.means2, 
        means=current.posterior.means2.1, 
        disps=current.posterior.disps2.1, 
        prior.location=current.posterior.mean.prior.location.1, 
        prior.scale=current.posterior.mean.prior.scale.1
      ) - 
      (log(current.posterior.means2.1) + 
         (log(current.posterior.means2.1) - log(proposed.posterior.means2.1))^2 / 
         (2*mean.proposal.scales2))
    
    current.posterior.means2.1[replace] <- proposed.posterior.means2.1[replace]
    accept.means2.1[replace] <- accept.means2.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.means2.2 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.means2.2), 
      sdlog=sqrt.mean.proposal.scales2
    )
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(
        n=samples2, 
        sample.means=sample.means2, 
        means=proposed.posterior.means2.2, 
        disps=current.posterior.disps2.2, 
        prior.location=current.posterior.mean.prior.location.2, 
        prior.scale=current.posterior.mean.prior.scale.2
      ) + 
      (log(proposed.posterior.means2.2) + 
         (log(proposed.posterior.means2.2) - log(current.posterior.means2.2))^2 / 
         (2*mean.proposal.scales2)) - 
      mean.conditional.log.posterior(
        n=samples2, 
        sample.means=sample.means2, 
        means=current.posterior.means2.2, 
        disps=current.posterior.disps2.2, 
        prior.location=current.posterior.mean.prior.location.2, 
        prior.scale=current.posterior.mean.prior.scale.2
      ) - 
      (log(current.posterior.means2.2) + 
         (log(current.posterior.means2.2) - log(proposed.posterior.means2.2))^2 / 
         (2*mean.proposal.scales2))
    
    current.posterior.means2.2[replace] <- proposed.posterior.means2.2[replace]
    accept.means2.2[replace] <- accept.means2.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.means2.3 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.means2.3), 
      sdlog=sqrt.mean.proposal.scales2
    )
    
    replace <- log(runif(genes)) <= 
      mean.conditional.log.posterior(
        n=samples2, 
        sample.means=sample.means2, 
        means=proposed.posterior.means2.3, 
        disps=current.posterior.disps2.3, 
        prior.location=current.posterior.mean.prior.location.3, 
        prior.scale=current.posterior.mean.prior.scale.3
      ) + 
      (log(proposed.posterior.means2.3) + 
         (log(proposed.posterior.means2.3) - log(current.posterior.means2.3))^2 / 
         (2*mean.proposal.scales2)) - 
      mean.conditional.log.posterior(
        n=samples2, 
        sample.means=sample.means2, 
        means=current.posterior.means2.3, 
        disps=current.posterior.disps2.3, 
        prior.location=current.posterior.mean.prior.location.3, 
        prior.scale=current.posterior.mean.prior.scale.3
      ) - 
      (log(current.posterior.means2.3) + 
         (log(current.posterior.means2.3) - log(proposed.posterior.means2.3))^2 / 
         (2*mean.proposal.scales2))
    
    current.posterior.means2.3[replace] <- proposed.posterior.means2.3[replace]
    accept.means2.3[replace] <- accept.means2.3[replace] + 1/chain.length
    
    # Metropolis updates for overall per-gene dispersions ####
    # Chain 1
    proposed.posterior.disps0.1 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.disps0.1), 
      sdlog=sqrt.disp.proposal.scales0
    )
    
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts, 
        n=samples0, 
        sample.means=sample.means0, 
        means=current.posterior.means0.1, 
        disps=proposed.posterior.disps0.1, 
        prior.location=current.posterior.disp.prior.location.1, 
        prior.scale=current.posterior.disp.prior.scale.1
      ) + 
      (log(proposed.posterior.disps0.1) + 
         (log(proposed.posterior.disps0.1) - log(current.posterior.disps0.1))^2 / 
         (2*disp.proposal.scales0)) - 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts, 
        n=samples0, 
        sample.means=sample.means0, 
        means=current.posterior.means0.1, 
        disps=current.posterior.disps0.1, 
        prior.location=current.posterior.disp.prior.location.1, 
        prior.scale=current.posterior.disp.prior.scale.1
      ) - 
      (log(current.posterior.disps0.1) + 
         (log(current.posterior.disps0.1) - log(proposed.posterior.disps0.1))^2 / 
         (2*disp.proposal.scales0))
    
    current.posterior.disps0.1[replace] <- proposed.posterior.disps0.1[replace]
    accept.disps0.1[replace] <- accept.disps0.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.disps0.2 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.disps0.2), 
      sdlog=sqrt.disp.proposal.scales0
    )
    
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts, 
        n=samples0, 
        sample.means=sample.means0, 
        means=current.posterior.means0.2, 
        disps=proposed.posterior.disps0.2, 
        prior.location=current.posterior.disp.prior.location.2, 
        prior.scale=current.posterior.disp.prior.scale.2
      ) + 
      (log(proposed.posterior.disps0.2) + 
         (log(proposed.posterior.disps0.2) - log(current.posterior.disps0.2))^2 / 
         (2*disp.proposal.scales0)) - 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts, 
        n=samples0, 
        sample.means=sample.means0, 
        means=current.posterior.means0.2, 
        disps=current.posterior.disps0.2, 
        prior.location=current.posterior.disp.prior.location.2, 
        prior.scale=current.posterior.disp.prior.scale.2
      ) - 
      (log(current.posterior.disps0.2) + 
         (log(current.posterior.disps0.2) - log(proposed.posterior.disps0.2))^2 / 
         (2*disp.proposal.scales0))
    
    current.posterior.disps0.2[replace] <- proposed.posterior.disps0.2[replace]
    accept.disps0.2[replace] <- accept.disps0.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.disps0.3 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.disps0.3), 
      sdlog=sqrt.disp.proposal.scales0
    )
    
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts, 
        n=samples0, 
        sample.means=sample.means0, 
        means=current.posterior.means0.3, 
        disps=proposed.posterior.disps0.3, 
        prior.location=current.posterior.disp.prior.location.3, 
        prior.scale=current.posterior.disp.prior.scale.3
      ) + 
      (log(proposed.posterior.disps0.3) + 
         (log(proposed.posterior.disps0.3) - log(current.posterior.disps0.3))^2 / 
         (2*disp.proposal.scales0)) - 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts, 
        n=samples0, 
        sample.means=sample.means0, 
        means=current.posterior.means0.3, 
        disps=current.posterior.disps0.3, 
        prior.location=current.posterior.disp.prior.location.3, 
        prior.scale=current.posterior.disp.prior.scale.3
      ) - 
      (log(current.posterior.disps0.3) + 
         (log(current.posterior.disps0.3) - log(proposed.posterior.disps0.3))^2 / 
         (2*disp.proposal.scales0))
    
    current.posterior.disps0.3[replace] <- proposed.posterior.disps0.3[replace]
    accept.disps0.3[replace] <- accept.disps0.3[replace] + 1/chain.length
    
    # Metropolis updates for group 1 per-gene dispersions ####
    # Chain 1
    proposed.posterior.disps1.1 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.disps1.1), 
      sdlog=sqrt.disp.proposal.scales1
    )
    
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts1, 
        n=samples1, 
        sample.means=sample.means1, 
        means=current.posterior.means1.1, 
        disps=proposed.posterior.disps1.1, 
        prior.location=current.posterior.disp.prior.location.1, 
        prior.scale=current.posterior.disp.prior.scale.1
      ) + 
      (log(proposed.posterior.disps1.1) + 
         (log(proposed.posterior.disps1.1) - log(current.posterior.disps1.1))^2 / 
         (2*disp.proposal.scales1)) - 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts1, 
        n=samples1, 
        sample.means=sample.means1, 
        means=current.posterior.means1.1, 
        disps=current.posterior.disps1.1, 
        prior.location=current.posterior.disp.prior.location.1, 
        prior.scale=current.posterior.disp.prior.scale.1
      ) - 
      (log(current.posterior.disps1.1) + 
         (log(current.posterior.disps1.1) - log(proposed.posterior.disps1.1))^2 / 
         (2*disp.proposal.scales1))
    
    current.posterior.disps1.1[replace] <- proposed.posterior.disps1.1[replace]
    accept.disps1.1[replace] <- accept.disps1.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.disps1.2 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.disps1.2), 
      sdlog=sqrt.disp.proposal.scales1
    )
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts1, 
        n=samples1, 
        sample.means=sample.means1, 
        means=current.posterior.means1.2, 
        disps=proposed.posterior.disps1.2, 
        prior.location=current.posterior.disp.prior.location.2, 
        prior.scale=current.posterior.disp.prior.scale.2
      ) + 
      (log(proposed.posterior.disps1.2) + 
         (log(proposed.posterior.disps1.2) - log(current.posterior.disps1.2))^2 / 
         (2*disp.proposal.scales1)) - 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts1, 
        n=samples1, 
        sample.means=sample.means1, 
        means=current.posterior.means1.2, 
        disps=current.posterior.disps1.2, 
        prior.location=current.posterior.disp.prior.location.2, 
        prior.scale=current.posterior.disp.prior.scale.2
      ) - 
      (log(current.posterior.disps1.2) + 
         (log(current.posterior.disps1.2) - log(proposed.posterior.disps1.2))^2 / 
         (2*disp.proposal.scales1))
    
    current.posterior.disps1.2[replace] <- proposed.posterior.disps1.2[replace]
    accept.disps1.2[replace] <- accept.disps1.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.disps1.3 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.disps1.3), 
      sdlog=sqrt.disp.proposal.scales1
    )
    
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts1, 
        n=samples1, 
        sample.means=sample.means1, 
        means=current.posterior.means1.3, 
        disps=proposed.posterior.disps1.3, 
        prior.location=current.posterior.disp.prior.location.3, 
        prior.scale=current.posterior.disp.prior.scale.3
      ) + 
      (log(proposed.posterior.disps1.3) + 
         (log(proposed.posterior.disps1.3) - log(current.posterior.disps1.3))^2 / 
         (2*disp.proposal.scales1)) - 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts1, 
        n=samples1, 
        sample.means=sample.means1, 
        means=current.posterior.means1.3, 
        disps=current.posterior.disps1.3, 
        prior.location=current.posterior.disp.prior.location.3, 
        prior.scale=current.posterior.disp.prior.scale.3
      ) - 
      (log(current.posterior.disps1.3) + 
         (log(current.posterior.disps1.3) - log(proposed.posterior.disps1.3))^2 / 
         (2*disp.proposal.scales1))
    
    current.posterior.disps1.3[replace] <- proposed.posterior.disps1.3[replace]
    accept.disps1.3[replace] <- accept.disps1.3[replace] + 1/chain.length
    
    # Metropolis updates for group 2 per-gene dispersions ####
    # Chain 1
    proposed.posterior.disps2.1 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.disps2.1), 
      sdlog=sqrt.disp.proposal.scales2
    )
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts2, 
        n=samples2, 
        sample.means=sample.means2, 
        means=current.posterior.means2.1, 
        disps=proposed.posterior.disps2.1, 
        prior.location=current.posterior.disp.prior.location.1, 
        prior.scale=current.posterior.disp.prior.scale.1
      ) + 
      (log(proposed.posterior.disps2.1) + 
         (log(proposed.posterior.disps2.1) - log(current.posterior.disps2.1))^2 / 
         (2*disp.proposal.scales2)) - 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts2, 
        n=samples2, 
        sample.means=sample.means2, 
        means=current.posterior.means2.1, 
        disps=current.posterior.disps2.1, 
        prior.location=current.posterior.disp.prior.location.1, 
        prior.scale=current.posterior.disp.prior.scale.1
      ) - 
      (log(current.posterior.disps2.1) + 
         (log(current.posterior.disps2.1) - log(proposed.posterior.disps2.1))^2 / 
         (2*disp.proposal.scales2))
    
    current.posterior.disps2.1[replace] <- proposed.posterior.disps2.1[replace]
    accept.disps2.1[replace] <- accept.disps2.1[replace] + 1/chain.length
    
    # Chain 2
    proposed.posterior.disps2.2 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.disps2.2), 
      sdlog=sqrt.disp.proposal.scales2
    )
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts2, 
        n=samples2, 
        sample.means=sample.means2, 
        means=current.posterior.means2.2, 
        disps=proposed.posterior.disps2.2, 
        prior.location=current.posterior.disp.prior.location.2, 
        prior.scale=current.posterior.disp.prior.scale.2
      ) + 
      (log(proposed.posterior.disps2.2) + 
         (log(proposed.posterior.disps2.2) - log(current.posterior.disps2.2))^2 / 
         (2*disp.proposal.scales2)) - 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts2, 
        n=samples2, 
        sample.means=sample.means2, 
        means=current.posterior.means2.2, 
        disps=current.posterior.disps2.2, 
        prior.location=current.posterior.disp.prior.location.2, 
        prior.scale=current.posterior.disp.prior.scale.2
      ) - 
      (log(current.posterior.disps2.2) + 
         (log(current.posterior.disps2.2) - log(proposed.posterior.disps2.2))^2 / 
         (2*disp.proposal.scales2))
    
    current.posterior.disps2.2[replace] <- proposed.posterior.disps2.2[replace]
    accept.disps2.2[replace] <- accept.disps2.2[replace] + 1/chain.length
    
    # Chain 3
    proposed.posterior.disps2.3 <- rlnorm(
      n=genes, 
      meanlog=log(current.posterior.disps2.3), 
      sdlog=sqrt.disp.proposal.scales2
    )
    
    replace <- log(runif(genes)) <= 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts2, 
        n=samples2, 
        sample.means=sample.means2, 
        means=current.posterior.means2.3, 
        disps=proposed.posterior.disps2.3, 
        prior.location=current.posterior.disp.prior.location.3, 
        prior.scale=current.posterior.disp.prior.scale.3
      ) + 
      (log(proposed.posterior.disps2.3) + 
         (log(proposed.posterior.disps2.3) - log(current.posterior.disps2.3))^2 / 
         (2*disp.proposal.scales2)) - 
      disp.conditional.log.posterior(
        genes=genes, 
        counts=counts2, 
        n=samples2, 
        sample.means=sample.means2, 
        means=current.posterior.means2.3, 
        disps=current.posterior.disps2.3, 
        prior.location=current.posterior.disp.prior.location.3, 
        prior.scale=current.posterior.disp.prior.scale.3
      ) - 
      (log(current.posterior.disps2.3) + 
         (log(current.posterior.disps2.3) - log(proposed.posterior.disps2.3))^2 / 
         (2*disp.proposal.scales2))
    
    current.posterior.disps2.3[replace] <- proposed.posterior.disps2.3[replace]
    accept.disps2.3[replace] <- accept.disps2.3[replace] + 1/chain.length
    
    # Gibbs update for prior location parameter for mean ####
    current.posterior.mean.prior.location.1 <- rnorm(
      n=1, 
      mean=prior.location.posterior.mean(
        hyperprior.mean=2.5, 
        hyperprior.var=20, 
        prior.scale=current.posterior.mean.prior.scale.1, 
        parameters0=current.posterior.means0.1, 
        parameters1=current.posterior.means1.1, 
        parameters2=current.posterior.means2.1, 
        z=current.posterior.indicators.1
      ), 
      sd=prior.location.posterior.sd(
        hyperprior.var=20, 
        prior.scale=current.posterior.mean.prior.scale.1, 
        z=current.posterior.indicators.1
      )
    )
    
    current.posterior.mean.prior.location.2 <- rnorm(
      n=1, 
      mean=prior.location.posterior.mean(
        hyperprior.mean=2.5, 
        hyperprior.var=20, 
        prior.scale=current.posterior.mean.prior.scale.2, 
        parameters0=current.posterior.means0.2, 
        parameters1=current.posterior.means1.2, 
        parameters2=current.posterior.means2.2, 
        z=current.posterior.indicators.2
      ), 
      sd=prior.location.posterior.sd(
        hyperprior.var=20, 
        prior.scale=current.posterior.mean.prior.scale.2, 
        z=current.posterior.indicators.2
      )
    )
    
    current.posterior.mean.prior.location.3 <- rnorm(
      n=1, 
      mean=prior.location.posterior.mean(
        hyperprior.mean=2.5, 
        hyperprior.var=20, 
        prior.scale=current.posterior.mean.prior.scale.3, 
        parameters0=current.posterior.means0.3, 
        parameters1=current.posterior.means1.3, 
        parameters2=current.posterior.means2.3, 
        z=current.posterior.indicators.3
      ), 
      sd=prior.location.posterior.sd(
        hyperprior.var=20, 
        prior.scale=current.posterior.mean.prior.scale.3, 
        z=current.posterior.indicators.3
      )
    )
    
    # Metropolis update for prior scale parameter for mean ####
    # Chain 1
    proposed.posterior.mean.prior.scale.1 <- rnorm(
      n=1, 
      mean=current.posterior.mean.prior.scale.1, 
      sd=mean.prior.scale.proposal.sd
    )
    if (log(runif(1)) <= 
        mean.prior.scale.log.posterior(
          prior.scale=proposed.posterior.mean.prior.scale.1, 
          prior.location=current.posterior.mean.prior.location.1, 
          g=genes, 
          means0=current.posterior.means0.1, 
          means1=current.posterior.means1.1, 
          means2=current.posterior.means2.1, 
          z=current.posterior.indicators.1
        ) - 
        mean.prior.scale.log.posterior(
          prior.scale=current.posterior.mean.prior.scale.1, 
          prior.location=current.posterior.mean.prior.location.1, 
          g=genes, 
          means0=current.posterior.means0.1, 
          means1=current.posterior.means1.1, 
          means2=current.posterior.means2.1, 
          z=current.posterior.indicators.1
        )) {
      current.posterior.mean.prior.scale.1 <- proposed.posterior.mean.prior.scale.1
      accept.mean.prior.scale.1 <- accept.mean.prior.scale.1 + 1/chain.length
    }
    
    # Chain 2
    proposed.posterior.mean.prior.scale.2 <- rnorm(
      n=1, 
      mean=current.posterior.mean.prior.scale.2, 
      sd=mean.prior.scale.proposal.sd
    )
    if (log(runif(1)) <= 
        mean.prior.scale.log.posterior(
          prior.scale=proposed.posterior.mean.prior.scale.2, 
          prior.location=current.posterior.mean.prior.location.2, 
          g=genes, 
          means0=current.posterior.means0.2, 
          means1=current.posterior.means1.2, 
          means2=current.posterior.means2.2, 
          z=current.posterior.indicators.2
        ) - 
        mean.prior.scale.log.posterior(
          prior.scale=current.posterior.mean.prior.scale.2, 
          prior.location=current.posterior.mean.prior.location.2, 
          g=genes, 
          means0=current.posterior.means0.2, 
          means1=current.posterior.means1.2, 
          means2=current.posterior.means2.2, 
          z=current.posterior.indicators.2
        )) {
      current.posterior.mean.prior.scale.2 <- proposed.posterior.mean.prior.scale.2
      accept.mean.prior.scale.2 <- accept.mean.prior.scale.2 + 1/chain.length
    }
    
    # Chain 3
    proposed.posterior.mean.prior.scale.3 <- rnorm(
      n=1, 
      mean=current.posterior.mean.prior.scale.3, 
      sd=mean.prior.scale.proposal.sd
    )
    if (log(runif(1)) <= 
        mean.prior.scale.log.posterior(
          prior.scale=proposed.posterior.mean.prior.scale.3, 
          prior.location=current.posterior.mean.prior.location.3, 
          g=genes, 
          means0=current.posterior.means0.3, 
          means1=current.posterior.means1.3, 
          means2=current.posterior.means2.3, 
          z=current.posterior.indicators.3
        ) - 
        mean.prior.scale.log.posterior(
          prior.scale=current.posterior.mean.prior.scale.3, 
          prior.location=current.posterior.mean.prior.location.3, 
          g=genes, 
          means0=current.posterior.means0.3, 
          means1=current.posterior.means1.3, 
          means2=current.posterior.means2.3, 
          z=current.posterior.indicators.3
        )) {
      current.posterior.mean.prior.scale.3 <- proposed.posterior.mean.prior.scale.3
      accept.mean.prior.scale.3 <- accept.mean.prior.scale.3 + 1/chain.length
    }
    
    # Gibbs update for prior location parameter for dispersion ####
    current.posterior.disp.prior.location.1 <- rnorm(
      n=1, 
      mean=prior.location.posterior.mean(
        hyperprior.mean=-2.5, 
        hyperprior.var=2, 
        prior.scale=current.posterior.disp.prior.scale.1, 
        parameters0=current.posterior.disps0.1, 
        parameters1=current.posterior.disps1.1, 
        parameters2=current.posterior.disps2.1, 
        z=current.posterior.indicators.1
      ), 
      sd=prior.location.posterior.sd(
        hyperprior.var=2, 
        prior.scale=current.posterior.disp.prior.scale.1, 
        z=current.posterior.indicators.1
      )
    )
    
    current.posterior.disp.prior.location.2 <- rnorm(
      n=1, 
      mean=prior.location.posterior.mean(
        hyperprior.mean=-2.5, 
        hyperprior.var=2, 
        prior.scale=current.posterior.disp.prior.scale.2, 
        parameters0=current.posterior.disps0.2, 
        parameters1=current.posterior.disps1.2, 
        parameters2=current.posterior.disps2.2, 
        z=current.posterior.indicators.2
      ), 
      sd=prior.location.posterior.sd(
        hyperprior.var=2, 
        prior.scale=current.posterior.disp.prior.scale.2, 
        z=current.posterior.indicators.2
      )
    )
    
    current.posterior.disp.prior.location.3 <- rnorm(
      n=1, 
      mean=prior.location.posterior.mean(
        hyperprior.mean=-2.5, 
        hyperprior.var=2, 
        prior.scale=current.posterior.disp.prior.scale.3, 
        parameters0=current.posterior.disps0.3, 
        parameters1=current.posterior.disps1.3, 
        parameters2=current.posterior.disps2.3, 
        z=current.posterior.indicators.3
      ), 
      sd=prior.location.posterior.sd(
        hyperprior.var=2, 
        prior.scale=current.posterior.disp.prior.scale.3, 
        z=current.posterior.indicators.3
      )
    )
    
    # Metropolis update for prior scale parameter for dispersion ####
    # Chain 1
    proposed.posterior.disp.prior.scale.1 <- rnorm(
      n=1, 
      mean=current.posterior.disp.prior.scale.1, 
      sd=disp.prior.scale.proposal.sd
    )
    if (log(runif(1)) <= 
        disp.prior.scale.log.posterior(
          prior.scale=proposed.posterior.disp.prior.scale.1, 
          prior.location=current.posterior.disp.prior.location.1, 
          g=genes, 
          disps0=current.posterior.disps0.1, 
          disps1=current.posterior.disps1.1, 
          disps2=current.posterior.disps2.1, 
          z=current.posterior.indicators.1
        ) - 
        disp.prior.scale.log.posterior(
          prior.scale=current.posterior.disp.prior.scale.1, 
          prior.location=current.posterior.disp.prior.location.1, 
          g=genes, 
          disps0=current.posterior.disps0.1, 
          disps1=current.posterior.disps1.1, 
          disps2=current.posterior.disps2.1, 
          z=current.posterior.indicators.1
        )) {
      current.posterior.disp.prior.scale.1 <- proposed.posterior.disp.prior.scale.1
      accept.disp.prior.scale.1 <- accept.disp.prior.scale.1 + 1/chain.length
    }
    
    # Chain 2
    proposed.posterior.disp.prior.scale.2 <- rnorm(
      n=1, 
      mean=current.posterior.disp.prior.scale.2, 
      sd=disp.prior.scale.proposal.sd
    )
    if (log(runif(1)) <= 
        disp.prior.scale.log.posterior(
          prior.scale=proposed.posterior.disp.prior.scale.2, 
          prior.location=current.posterior.disp.prior.location.2, 
          g=genes, 
          disps0=current.posterior.disps0.2, 
          disps1=current.posterior.disps1.2, 
          disps2=current.posterior.disps2.2, 
          z=current.posterior.indicators.2
        ) - 
        disp.prior.scale.log.posterior(
          prior.scale=current.posterior.disp.prior.scale.2, 
          prior.location=current.posterior.disp.prior.location.2, 
          g=genes, 
          disps0=current.posterior.disps0.2, 
          disps1=current.posterior.disps1.2, 
          disps2=current.posterior.disps2.2, 
          z=current.posterior.indicators.2
        )) {
      current.posterior.disp.prior.scale.2 <- proposed.posterior.disp.prior.scale.2
      accept.disp.prior.scale.2 <- accept.disp.prior.scale.2 + 1/chain.length
    }
    
    # Chain 3
    proposed.posterior.disp.prior.scale.3 <- rnorm(
      n=1, 
      mean=current.posterior.disp.prior.scale.3, 
      sd=disp.prior.scale.proposal.sd
    )
    if (log(runif(1)) <= 
        disp.prior.scale.log.posterior(
          prior.scale=proposed.posterior.disp.prior.scale.3, 
          prior.location=current.posterior.disp.prior.location.3, 
          g=genes, 
          disps0=current.posterior.disps0.3, 
          disps1=current.posterior.disps1.3, 
          disps2=current.posterior.disps2.3, 
          z=current.posterior.indicators.3
        ) - 
        disp.prior.scale.log.posterior(
          prior.scale=current.posterior.disp.prior.scale.3, 
          prior.location=current.posterior.disp.prior.location.3, 
          g=genes, 
          disps0=current.posterior.disps0.3, 
          disps1=current.posterior.disps1.3, 
          disps2=current.posterior.disps2.3, 
          z=current.posterior.indicators.3
        )) {
      current.posterior.disp.prior.scale.3 <- proposed.posterior.disp.prior.scale.3
      accept.disp.prior.scale.3 <- accept.disp.prior.scale.3 + 1/chain.length
    }
    
    # Gibbs updates for per-gene mixture components ####
    current.posterior.indicators.1 <- rbinom(
      n=genes, 
      size=1, 
      prob=posterior.indicator.probabilities(
        genes=genes, 
        counts=counts, 
        counts1=counts1, 
        counts2=counts2, 
        n=samples0, 
        n1=samples1, 
        n2=samples2, 
        sample.means0=sample.means0, 
        sample.means1=sample.means1, 
        sample.means2=sample.means2,
        means0=current.posterior.means0.1, 
        means1=current.posterior.means1.1, 
        means2=current.posterior.means2.1, 
        disps0=current.posterior.disps0.1, 
        disps1=current.posterior.disps1.1, 
        disps2=current.posterior.disps2.1, 
        mean.prior.location=current.posterior.mean.prior.location.1, 
        mean.prior.scale=current.posterior.mean.prior.scale.1, 
        disp.prior.location=current.posterior.disp.prior.location.1, 
        disp.prior.scale=current.posterior.disp.prior.scale.1, 
        lambda=current.posterior.proportion.1
      )
    )
    
    current.posterior.indicators.2 <- rbinom(
      n=genes, 
      size=1, 
      prob=posterior.indicator.probabilities(
        genes=genes, 
        counts=counts, 
        counts1=counts1, 
        counts2=counts2, 
        n=samples0, 
        n1=samples1, 
        n2=samples2, 
        sample.means0=sample.means0, 
        sample.means1=sample.means1, 
        sample.means2=sample.means2,
        means0=current.posterior.means0.2, 
        means1=current.posterior.means1.2, 
        means2=current.posterior.means2.2, 
        disps0=current.posterior.disps0.2, 
        disps1=current.posterior.disps1.2, 
        disps2=current.posterior.disps2.2, 
        mean.prior.location=current.posterior.mean.prior.location.2, 
        mean.prior.scale=current.posterior.mean.prior.scale.2, 
        disp.prior.location=current.posterior.disp.prior.location.2, 
        disp.prior.scale=current.posterior.disp.prior.scale.2, 
        lambda=current.posterior.proportion.2
      )
    )
    
    current.posterior.indicators.3 <- rbinom(
      n=genes, 
      size=1, 
      prob=posterior.indicator.probabilities(
        genes=genes, 
        counts=counts, 
        counts1=counts1, 
        counts2=counts2, 
        n=samples0, 
        n1=samples1, 
        n2=samples2, 
        sample.means0=sample.means0, 
        sample.means1=sample.means1, 
        sample.means2=sample.means2,
        means0=current.posterior.means0.3, 
        means1=current.posterior.means1.3, 
        means2=current.posterior.means2.3, 
        disps0=current.posterior.disps0.3, 
        disps1=current.posterior.disps1.3, 
        disps2=current.posterior.disps2.3, 
        mean.prior.location=current.posterior.mean.prior.location.3, 
        mean.prior.scale=current.posterior.mean.prior.scale.3, 
        disp.prior.location=current.posterior.disp.prior.location.3, 
        disp.prior.scale=current.posterior.disp.prior.scale.3, 
        lambda=current.posterior.proportion.3
      )
    )
    
    # Gibbs update for mixture proportion ####
    current.posterior.proportion.1 <- rbeta(
      n=1, 
      shape1=1 + sum(current.posterior.indicators.1), 
      shape2=1 + genes - sum(current.posterior.indicators.1)
    )
    
    current.posterior.proportion.2 <- rbeta(
      n=1, 
      shape1=1 + sum(current.posterior.indicators.2), 
      shape2=1 + genes - sum(current.posterior.indicators.2)
    )
    
    current.posterior.proportion.3 <- rbeta(
      n=1, 
      shape1=1 + sum(current.posterior.indicators.3), 
      shape2=1 + genes - sum(current.posterior.indicators.3)
    )
    
    # Update posterior samples ####
    if (iter/thin==round(iter/thin)) {
      posterior.means0.1[iter/thin,] <- current.posterior.means0.1
      posterior.means1.1[iter/thin,] <- current.posterior.means1.1
      posterior.means2.1[iter/thin,] <- current.posterior.means2.1
      posterior.means0.2[iter/thin,] <- current.posterior.means0.2
      posterior.means1.2[iter/thin,] <- current.posterior.means1.2
      posterior.means2.2[iter/thin,] <- current.posterior.means2.2
      posterior.means0.3[iter/thin,] <- current.posterior.means0.3
      posterior.means1.3[iter/thin,] <- current.posterior.means1.3
      posterior.means2.3[iter/thin,] <- current.posterior.means2.3
      posterior.disps0.1[iter/thin,] <- current.posterior.disps0.1
      posterior.disps1.1[iter/thin,] <- current.posterior.disps1.1
      posterior.disps2.1[iter/thin,] <- current.posterior.disps2.1
      posterior.disps0.2[iter/thin,] <- current.posterior.disps0.2
      posterior.disps1.2[iter/thin,] <- current.posterior.disps1.2
      posterior.disps2.2[iter/thin,] <- current.posterior.disps2.2
      posterior.disps0.3[iter/thin,] <- current.posterior.disps0.3
      posterior.disps1.3[iter/thin,] <- current.posterior.disps1.3
      posterior.disps2.3[iter/thin,] <- current.posterior.disps2.3
      posterior.mean.prior.location.1[iter/thin] <- 
        current.posterior.mean.prior.location.1
      posterior.mean.prior.location.2[iter/thin] <- 
        current.posterior.mean.prior.location.2
      posterior.mean.prior.location.3[iter/thin] <- 
        current.posterior.mean.prior.location.3
      posterior.mean.prior.scale.1[iter/thin] <- 
        current.posterior.mean.prior.scale.1
      posterior.mean.prior.scale.2[iter/thin] <- 
        current.posterior.mean.prior.scale.2
      posterior.mean.prior.scale.3[iter/thin] <- 
        current.posterior.mean.prior.scale.3
      posterior.disp.prior.location.1[iter/thin] <- 
        current.posterior.disp.prior.location.1
      posterior.disp.prior.location.2[iter/thin] <- 
        current.posterior.disp.prior.location.2
      posterior.disp.prior.location.3[iter/thin] <- 
        current.posterior.disp.prior.location.3
      posterior.disp.prior.scale.1[iter/thin] <- current.posterior.disp.prior.scale.1
      posterior.disp.prior.scale.2[iter/thin] <- current.posterior.disp.prior.scale.2
      posterior.disp.prior.scale.3[iter/thin] <- current.posterior.disp.prior.scale.3
      posterior.indicators.1[iter/thin,] <- current.posterior.indicators.1
      posterior.indicators.2[iter/thin,] <- current.posterior.indicators.2
      posterior.indicators.3[iter/thin,] <- current.posterior.indicators.3
      posterior.proportion.1[iter/thin] <- current.posterior.proportion.1
      posterior.proportion.2[iter/thin] <- current.posterior.proportion.2
      posterior.proportion.3[iter/thin] <- current.posterior.proportion.3
    }
    
  }
  
  return(list("chain.length"=chain.length, 
              "thin"=thin, 
              "inits1"=inits1, 
              "inits2"=inits2, 
              "inits3"=inits3, 
              "mean.proposal.scales0"=mean.proposal.scales0, 
              "mean.proposal.scales1"=mean.proposal.scales1, 
              "mean.proposal.scales2"=mean.proposal.scales2, 
              "disp.proposal.scales0"=disp.proposal.scales0, 
              "disp.proposal.scales1"=disp.proposal.scales1, 
              "disp.proposal.scales2"=disp.proposal.scales2, 
              "mean.prior.scale.proposal.sd"=mean.prior.scale.proposal.sd, 
              "disp.prior.scale.proposal.sd"=disp.prior.scale.proposal.sd, 
              "accept.means0.1"=accept.means0.1, 
              "accept.means1.1"=accept.means1.1, 
              "accept.means2.1"=accept.means2.1, 
              "accept.means0.2"=accept.means0.2, 
              "accept.means1.2"=accept.means1.2, 
              "accept.means2.2"=accept.means2.2, 
              "accept.means0.3"=accept.means0.3, 
              "accept.means1.3"=accept.means1.3, 
              "accept.means2.3"=accept.means2.3, 
              "accept.disps0.1"=accept.disps0.1, 
              "accept.disps1.1"=accept.disps1.1, 
              "accept.disps2.1"=accept.disps2.1, 
              "accept.disps0.2"=accept.disps0.2, 
              "accept.disps1.2"=accept.disps1.2, 
              "accept.disps2.2"=accept.disps2.2, 
              "accept.disps0.3"=accept.disps0.3, 
              "accept.disps1.3"=accept.disps1.3, 
              "accept.disps2.3"=accept.disps2.3, 
              "accept.mean.prior.scale.1"=accept.mean.prior.scale.1, 
              "accept.mean.prior.scale.2"=accept.mean.prior.scale.2, 
              "accept.mean.prior.scale.3"=accept.mean.prior.scale.3, 
              "accept.disp.prior.scale.1"=accept.disp.prior.scale.1, 
              "accept.disp.prior.scale.2"=accept.disp.prior.scale.2, 
              "accept.disp.prior.scale.3"=accept.disp.prior.scale.3, 
              "posterior.means0.1"=posterior.means0.1, 
              "posterior.means1.1"=posterior.means1.1, 
              "posterior.means2.1"=posterior.means2.1, 
              "posterior.means0.2"=posterior.means0.2, 
              "posterior.means1.2"=posterior.means1.2, 
              "posterior.means2.2"=posterior.means2.2, 
              "posterior.means0.3"=posterior.means0.3, 
              "posterior.means1.3"=posterior.means1.3, 
              "posterior.means2.3"=posterior.means2.3, 
              "posterior.disps0.1"=posterior.disps0.1, 
              "posterior.disps1.1"=posterior.disps1.1, 
              "posterior.disps2.1"=posterior.disps2.1, 
              "posterior.disps0.2"=posterior.disps0.2, 
              "posterior.disps1.2"=posterior.disps1.2, 
              "posterior.disps2.2"=posterior.disps2.2, 
              "posterior.disps0.3"=posterior.disps0.3, 
              "posterior.disps1.3"=posterior.disps1.3, 
              "posterior.disps2.3"=posterior.disps2.3, 
              "posterior.mean.prior.location.1"=posterior.mean.prior.location.1, 
              "posterior.mean.prior.location.2"=posterior.mean.prior.location.2, 
              "posterior.mean.prior.location.3"=posterior.mean.prior.location.3, 
              "posterior.mean.prior.scale.1"=posterior.mean.prior.scale.1, 
              "posterior.mean.prior.scale.2"=posterior.mean.prior.scale.2, 
              "posterior.mean.prior.scale.3"=posterior.mean.prior.scale.3, 
              "posterior.disp.prior.location.1"=posterior.disp.prior.location.1, 
              "posterior.disp.prior.location.2"=posterior.disp.prior.location.2, 
              "posterior.disp.prior.location.3"=posterior.disp.prior.location.3, 
              "posterior.disp.prior.scale.1"=posterior.disp.prior.scale.1, 
              "posterior.disp.prior.scale.2"=posterior.disp.prior.scale.2, 
              "posterior.disp.prior.scale.3"=posterior.disp.prior.scale.3, 
              "posterior.indicators.1"=posterior.indicators.1, 
              "posterior.indicators.2"=posterior.indicators.2, 
              "posterior.indicators.3"=posterior.indicators.3, 
              "posterior.proportion.1"=posterior.proportion.1, 
              "posterior.proportion.2"=posterior.proportion.2, 
              "posterior.proportion.3"=posterior.proportion.3))
  
}
