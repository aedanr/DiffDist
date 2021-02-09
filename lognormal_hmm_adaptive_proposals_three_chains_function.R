ln_hmm_adapt_3_chains <- function(counts, 
                                  groups, 
                                  chain.length=2000, 
                                  initial.chain.length=100, 
                                  adapt.chain.length=100) {
  
  require(coda)
  genes <- ncol(counts)
  counts1 <- counts[groups==1,]
  counts2 <- counts[groups==2,]
  samples0 <- nrow(counts)
  samples1 <- nrow(counts1)
  samples2 <- nrow(counts2)
  sample.means0 <- pmax(colMeans(counts), 0.01)
  sample.means1 <- pmax(colMeans(counts1), 0.01)
  sample.means2 <- pmax(colMeans(counts2), 0.01)
  sample.vars0 <- apply(counts, 2, var)
  sample.vars1 <- apply(counts1, 2, var)
  sample.vars2 <- apply(counts2, 2, var)
  rm(counts1, counts2)
  sample.disps0 <- pmax(sample.vars0-sample.means0, 0.01) / pmax(sample.means0^2, 0.1)
  sample.disps1 <- pmax(sample.vars1-sample.means1, 0.01) / pmax(sample.means1^2, 0.1)
  sample.disps2 <- pmax(sample.vars2-sample.means2, 0.01) / pmax(sample.means2^2, 0.1)
  rm(sample.vars0, sample.vars1, sample.vars2)
  
  # Set proposal scales
  message("initialising...")
  adapt.mcmc <- ln_hmm_1_chain(
    counts=counts, 
    groups=groups, 
    chain.length=initial.chain.length, 
    thin=1, 
    inits=list("means0" = sample.means0, 
               "means1" = sample.means1, 
               "means2" = sample.means2, 
               "disps0" = sample.disps0, 
               "disps1" = sample.disps1, 
               "disps2" = sample.disps2, 
               "mean.prior.location" = 1, 
               "disp.prior.location" = 1, 
               "mean.prior.scale" = 1, 
               "disp.prior.scale" = 1), 
    mean.proposal.scales0=rep(0.2, ncol(counts)), 
    mean.proposal.scales1=rep(0.2, ncol(counts)), 
    mean.proposal.scales2=rep(0.2, ncol(counts)), 
    disp.proposal.scales0=rep(0.5, ncol(counts)), 
    disp.proposal.scales1=rep(0.5, ncol(counts)), 
    disp.proposal.scales2=rep(0.5, ncol(counts)), 
    mean.prior.scale.proposal.sd=0.1, 
    disp.prior.scale.proposal.sd=0.4
  )
  rm(sample.means0, sample.means1, sample.means2, 
     sample.disps0, sample.disps1, sample.disps2)
  
  inits1.adapt <- list(
    "means0"=apply(adapt.mcmc$posterior.means0, 2, min)/5, 
    "means1"=apply(adapt.mcmc$posterior.means1, 2, min)/5, 
    "means2"=apply(adapt.mcmc$posterior.means2, 2, min)/5, 
    "disps0"=apply(adapt.mcmc$posterior.disps0, 2, min)/10, 
    "disps1"=apply(adapt.mcmc$posterior.disps1, 2, min)/10, 
    "disps2"=apply(adapt.mcmc$posterior.disps2, 2, min)/10, 
    "mean.prior.location"=min(adapt.mcmc$posterior.mean.prior.location)/2,
    "mean.prior.scale"=min(adapt.mcmc$posterior.mean.prior.scale)/2,
    "disp.prior.location"=min(adapt.mcmc$posterior.disp.prior.location)/2, 
    "disp.prior.scale"=min(adapt.mcmc$posterior.disp.prior.scale)/2
  )
  inits2.adapt <- list(
    "means0"=colMeans(adapt.mcmc$posterior.means0), 
    "means1"=colMeans(adapt.mcmc$posterior.means1), 
    "means2"=colMeans(adapt.mcmc$posterior.means2), 
    "disps0"=colMeans(adapt.mcmc$posterior.disps0), 
    "disps1"=colMeans(adapt.mcmc$posterior.disps1), 
    "disps2"=colMeans(adapt.mcmc$posterior.disps2), 
    "mean.prior.location"=mean(adapt.mcmc$posterior.mean.prior.location),
    "mean.prior.scale"=mean(adapt.mcmc$posterior.mean.prior.scale),
    "disp.prior.location"=mean(adapt.mcmc$posterior.disp.prior.location), 
    "disp.prior.scale"=mean(adapt.mcmc$posterior.disp.prior.scale)
  )
  inits3.adapt <- list(
    "means0"=apply(adapt.mcmc$posterior.means0, 2, max)*5, 
    "means1"=apply(adapt.mcmc$posterior.means1, 2, max)*5, 
    "means2"=apply(adapt.mcmc$posterior.means2, 2, max)*5, 
    "disps0"=apply(adapt.mcmc$posterior.disps0, 2, max)*10, 
    "disps1"=apply(adapt.mcmc$posterior.disps1, 2, max)*10, 
    "disps2"=apply(adapt.mcmc$posterior.disps2, 2, max)*10, 
    "mean.prior.location"=max(adapt.mcmc$posterior.mean.prior.location)*2,
    "mean.prior.scale"=max(adapt.mcmc$posterior.mean.prior.scale)*2,
    "disp.prior.location"=max(adapt.mcmc$posterior.disp.prior.location)*2, 
    "disp.prior.scale"=max(adapt.mcmc$posterior.disp.prior.scale)*2
  )
  
  message("optimising proposal distributions...")
  current.mean.proposal.scales0 <- 
    adapt.mcmc$mean.proposal.scales0 * (1 + 2.25 * (adapt.mcmc$accept.means0 - 0.44))
  current.mean.proposal.scales1 <- 
    adapt.mcmc$mean.proposal.scales1 * (1 + 2.25 * (adapt.mcmc$accept.means1 - 0.44))
  current.mean.proposal.scales2 <- 
    adapt.mcmc$mean.proposal.scales2 * (1 + 2.25 * (adapt.mcmc$accept.means2 - 0.44))
  current.disp.proposal.scales0 <- 
    adapt.mcmc$disp.proposal.scales0 * (1 + 2.25 * (adapt.mcmc$accept.disps0 - 0.44))
  current.disp.proposal.scales1 <- 
    adapt.mcmc$disp.proposal.scales1 * (1 + 2.25 * (adapt.mcmc$accept.disps1 - 0.44))
  current.disp.proposal.scales2 <- 
    adapt.mcmc$disp.proposal.scales2 * (1 + 2.25 * (adapt.mcmc$accept.disps2 - 0.44))
  current.mean.prior.scale.proposal.sd <- 
    adapt.mcmc$mean.prior.scale.proposal.sd * 
    (1 + 2.25 * (adapt.mcmc$accept.mean.prior.scale - 0.44))
  current.disp.prior.scale.proposal.sd <- 
    adapt.mcmc$disp.prior.scale.proposal.sd * 
    (1 + 2.25 * (adapt.mcmc$accept.disp.prior.scale - 0.44))
  
  adapt.mcmc <- ln_hmm_3_chains(
    counts=counts, 
    groups=groups, 
    chain.length=adapt.chain.length, 
    thin=adapt.chain.length, 
    inits1=inits1.adapt, 
    inits2=inits2.adapt, 
    inits3=inits3.adapt, 
    mean.proposal.scales0=current.mean.proposal.scales0, 
    mean.proposal.scales1=current.mean.proposal.scales1, 
    mean.proposal.scales2=current.mean.proposal.scales2, 
    disp.proposal.scales0=current.disp.proposal.scales0, 
    disp.proposal.scales1=current.disp.proposal.scales1, 
    disp.proposal.scales2=current.disp.proposal.scales2, 
    mean.prior.scale.proposal.sd=current.mean.prior.scale.proposal.sd, 
    disp.prior.scale.proposal.sd=current.disp.prior.scale.proposal.sd
  )
  rm(inits1.adapt, inits2.adapt, inits3.adapt)
  
  current.means0.1 <- adapt.mcmc$posterior.means0.1
  current.means1.1 <- adapt.mcmc$posterior.means1.1
  current.means2.1 <- adapt.mcmc$posterior.means2.1
  current.means0.2 <- adapt.mcmc$posterior.means0.2
  current.means1.2 <- adapt.mcmc$posterior.means1.2
  current.means2.2 <- adapt.mcmc$posterior.means2.2
  current.means0.3 <- adapt.mcmc$posterior.means0.3
  current.means1.3 <- adapt.mcmc$posterior.means1.3
  current.means2.3 <- adapt.mcmc$posterior.means2.3
  current.disps0.1 <- adapt.mcmc$posterior.disps0.1
  current.disps1.1 <- adapt.mcmc$posterior.disps1.1
  current.disps2.1 <- adapt.mcmc$posterior.disps2.1
  current.disps0.2 <- adapt.mcmc$posterior.disps0.2
  current.disps1.2 <- adapt.mcmc$posterior.disps1.2
  current.disps2.2 <- adapt.mcmc$posterior.disps2.2
  current.disps0.3 <- adapt.mcmc$posterior.disps0.3
  current.disps1.3 <- adapt.mcmc$posterior.disps1.3
  current.disps2.3 <- adapt.mcmc$posterior.disps2.3
  current.mean.prior.location.1 <- adapt.mcmc$posterior.mean.prior.location.1
  current.mean.prior.location.2 <- adapt.mcmc$posterior.mean.prior.location.2
  current.mean.prior.location.3 <- adapt.mcmc$posterior.mean.prior.location.3
  current.mean.prior.scale.1 <- adapt.mcmc$posterior.mean.prior.scale.1
  current.mean.prior.scale.2 <- adapt.mcmc$posterior.mean.prior.scale.2
  current.mean.prior.scale.3 <- adapt.mcmc$posterior.mean.prior.scale.3
  current.disp.prior.location.1 <- adapt.mcmc$posterior.disp.prior.location.1
  current.disp.prior.location.2 <- adapt.mcmc$posterior.disp.prior.location.2
  current.disp.prior.location.3 <- adapt.mcmc$posterior.disp.prior.location.3
  current.disp.prior.scale.1 <- adapt.mcmc$posterior.disp.prior.scale.1
  current.disp.prior.scale.2 <- adapt.mcmc$posterior.disp.prior.scale.2
  current.disp.prior.scale.3 <- adapt.mcmc$posterior.disp.prior.scale.3
  
  accept.means0 <- 
    (adapt.mcmc$accept.means0.1 + adapt.mcmc$accept.means0.2 + 
       adapt.mcmc$accept.means0.3)/3
  accept.means1 <- 
    (adapt.mcmc$accept.means1.1 + adapt.mcmc$accept.means1.2 + 
       adapt.mcmc$accept.means1.3)/3
  accept.means2 <- 
    (adapt.mcmc$accept.means2.1 + adapt.mcmc$accept.means2.2 + 
       adapt.mcmc$accept.means2.3)/3
  accept.disps0 <- 
    (adapt.mcmc$accept.disps0.1 + adapt.mcmc$accept.disps0.2 + 
       adapt.mcmc$accept.disps0.3)/3
  accept.disps1 <- 
    (adapt.mcmc$accept.disps1.1 + adapt.mcmc$accept.disps1.2 + 
       adapt.mcmc$accept.disps1.3)/3
  accept.disps2 <- 
    (adapt.mcmc$accept.disps2.1 + adapt.mcmc$accept.disps2.2 + 
       adapt.mcmc$accept.disps2.3)/3
  accept.mean.prior.scale <- mean(c(adapt.mcmc$accept.mean.prior.scale.1, 
                                    adapt.mcmc$accept.mean.prior.scale.2, 
                                    adapt.mcmc$accept.mean.prior.scale.3))
  accept.disp.prior.scale <- mean(c(adapt.mcmc$accept.disp.prior.scale.1, 
                                    adapt.mcmc$accept.disp.prior.scale.2, 
                                    adapt.mcmc$accept.disp.prior.scale.3))
  adapt.runs <- 1
  
  while (sum(accept.means0<0.24 | accept.means0>0.64 | 
             accept.means1<0.24 | accept.means1>0.64 | 
             accept.means2<0.24 | accept.means2>0.64 | 
             accept.disps0<0.24 | accept.disps0>0.64 | 
             accept.disps1<0.24 | accept.disps1>0.64 | 
             accept.disps2<0.24 | accept.disps2>0.64 | 
             accept.mean.prior.scale<0.24 | accept.mean.prior.scale>0.64 | 
             accept.disp.prior.scale<0.24 | accept.disp.prior.scale>0.64) > 0) {
    
    current.mean.proposal.scales0 <- 
      current.mean.proposal.scales0 * (1 + 2.25 * (accept.means0-0.44))
    current.mean.proposal.scales1 <- 
      current.mean.proposal.scales1 * (1 + 2.25 * (accept.means1-0.44))
    current.mean.proposal.scales2 <- 
      current.mean.proposal.scales2 * (1 + 2.25 * (accept.means2-0.44))
    current.disp.proposal.scales0 <- 
      current.disp.proposal.scales0 * (1 + 2.25 * (accept.disps0-0.44))
    current.disp.proposal.scales1 <- 
      current.disp.proposal.scales1 * (1 + 2.25 * (accept.disps1-0.44))
    current.disp.proposal.scales2 <- 
      current.disp.proposal.scales2 * (1 + 2.25 * (accept.disps2-0.44))
    current.mean.prior.scale.proposal.sd <- 
      current.mean.prior.scale.proposal.sd*(1+2.25*(accept.mean.prior.scale-0.44))
    current.disp.prior.scale.proposal.sd <- 
      current.disp.prior.scale.proposal.sd*(1+2.25*(accept.disp.prior.scale-0.44))
    
    adapt.mcmc <- ln_hmm_3_chains(
      counts=counts, 
      groups=groups, 
      chain.length=adapt.chain.length, 
      thin=adapt.chain.length, 
      inits1=list("means0"=current.means0.1, 
                  "means1"=current.means1.1, 
                  "means2"=current.means2.1, 
                  "disps0"=current.disps0.1, 
                  "disps1"=current.disps1.1, 
                  "disps2"=current.disps2.1, 
                  "mean.prior.location"=current.mean.prior.location.1, 
                  "mean.prior.scale"=current.mean.prior.scale.1, 
                  "disp.prior.location"=current.disp.prior.location.1, 
                  "disp.prior.scale"=current.disp.prior.scale.1), 
      inits2=list("means0"=current.means0.2, 
                  "means1"=current.means1.2, 
                  "means2"=current.means2.2, 
                  "disps0"=current.disps0.2, 
                  "disps1"=current.disps1.2, 
                  "disps2"=current.disps2.2, 
                  "mean.prior.location"=current.mean.prior.location.2, 
                  "mean.prior.scale"=current.mean.prior.scale.2, 
                  "disp.prior.location"=current.disp.prior.location.2, 
                  "disp.prior.scale"=current.disp.prior.scale.2), 
      inits3=list("means0"=current.means0.3, 
                  "means1"=current.means1.3, 
                  "means2"=current.means2.3, 
                  "disps0"=current.disps0.3, 
                  "disps1"=current.disps1.3, 
                  "disps2"=current.disps2.3, 
                  "mean.prior.location"=current.mean.prior.location.3, 
                  "mean.prior.scale"=current.mean.prior.scale.3, 
                  "disp.prior.location"=current.disp.prior.location.3, 
                  "disp.prior.scale"=current.disp.prior.scale.3), 
      mean.proposal.scales0=current.mean.proposal.scales0, 
      mean.proposal.scales1=current.mean.proposal.scales1, 
      mean.proposal.scales2=current.mean.proposal.scales2, 
      disp.proposal.scales0=current.disp.proposal.scales0, 
      disp.proposal.scales1=current.disp.proposal.scales1, 
      disp.proposal.scales2=current.disp.proposal.scales2, 
      mean.prior.scale.proposal.sd=current.mean.prior.scale.proposal.sd, 
      disp.prior.scale.proposal.sd=current.disp.prior.scale.proposal.sd
    )
    
    current.means0.1 <- adapt.mcmc$posterior.means0.1
    current.means1.1 <- adapt.mcmc$posterior.means1.1
    current.means2.1 <- adapt.mcmc$posterior.means2.1
    current.means0.2 <- adapt.mcmc$posterior.means0.2
    current.means1.2 <- adapt.mcmc$posterior.means1.2
    current.means2.2 <- adapt.mcmc$posterior.means2.2
    current.means0.3 <- adapt.mcmc$posterior.means0.3
    current.means1.3 <- adapt.mcmc$posterior.means1.3
    current.means2.3 <- adapt.mcmc$posterior.means2.3
    current.disps0.1 <- adapt.mcmc$posterior.disps0.1
    current.disps1.1 <- adapt.mcmc$posterior.disps1.1
    current.disps2.1 <- adapt.mcmc$posterior.disps2.1
    current.disps0.2 <- adapt.mcmc$posterior.disps0.2
    current.disps1.2 <- adapt.mcmc$posterior.disps1.2
    current.disps2.2 <- adapt.mcmc$posterior.disps2.2
    current.disps0.3 <- adapt.mcmc$posterior.disps0.3
    current.disps1.3 <- adapt.mcmc$posterior.disps1.3
    current.disps2.3 <- adapt.mcmc$posterior.disps2.3
    current.mean.prior.location.1 <- adapt.mcmc$posterior.mean.prior.location.1
    current.mean.prior.location.2 <- adapt.mcmc$posterior.mean.prior.location.2
    current.mean.prior.location.3 <- adapt.mcmc$posterior.mean.prior.location.3
    current.mean.prior.scale.1 <- adapt.mcmc$posterior.mean.prior.scale.1
    current.mean.prior.scale.2 <- adapt.mcmc$posterior.mean.prior.scale.2
    current.mean.prior.scale.3 <- adapt.mcmc$posterior.mean.prior.scale.3
    current.disp.prior.location.1 <- adapt.mcmc$posterior.disp.prior.location.1
    current.disp.prior.location.2 <- adapt.mcmc$posterior.disp.prior.location.2
    current.disp.prior.location.3 <- adapt.mcmc$posterior.disp.prior.location.3
    current.disp.prior.scale.1 <- adapt.mcmc$posterior.disp.prior.scale.1
    current.disp.prior.scale.2 <- adapt.mcmc$posterior.disp.prior.scale.2
    current.disp.prior.scale.3 <- adapt.mcmc$posterior.disp.prior.scale.3
    
    accept.means0 <- 
      (adapt.mcmc$accept.means0.1 + adapt.mcmc$accept.means0.2 + 
         adapt.mcmc$accept.means0.3)/3
    accept.means1 <- 
      (adapt.mcmc$accept.means1.1 + adapt.mcmc$accept.means1.2 + 
         adapt.mcmc$accept.means1.3)/3
    accept.means2 <- 
      (adapt.mcmc$accept.means2.1 + adapt.mcmc$accept.means2.2 + 
         adapt.mcmc$accept.means2.3)/3
    accept.disps0 <- 
      (adapt.mcmc$accept.disps0.1 + adapt.mcmc$accept.disps0.2 + 
         adapt.mcmc$accept.disps0.3)/3
    accept.disps1 <- 
      (adapt.mcmc$accept.disps1.1 + adapt.mcmc$accept.disps1.2 + 
         adapt.mcmc$accept.disps1.3)/3
    accept.disps2 <- 
      (adapt.mcmc$accept.disps2.1 + adapt.mcmc$accept.disps2.2 + 
         adapt.mcmc$accept.disps2.3)/3
    accept.mean.prior.scale <- mean(c(adapt.mcmc$accept.mean.prior.scale.1, 
                                      adapt.mcmc$accept.mean.prior.scale.2, 
                                      adapt.mcmc$accept.mean.prior.scale.3))
    accept.disp.prior.scale <- mean(c(adapt.mcmc$accept.disp.prior.scale.1, 
                                      adapt.mcmc$accept.disp.prior.scale.2, 
                                      adapt.mcmc$accept.disp.prior.scale.3))
    
    adapt.runs <- adapt.runs+1
  }
  rm(adapt.mcmc)
  
  inits1 <- list("means0"=current.means0.1, "means1"=current.means1.1, 
                 "means2"=current.means2.1, "disps0"=current.disps0.1, 
                 "disps1"=current.disps1.1, "disps2"=current.disps2.1, 
                 "mean.prior.location"=current.mean.prior.location.1, 
                 "mean.prior.scale"=current.mean.prior.scale.1, 
                 "disp.prior.location"=current.disp.prior.location.1, 
                 "disp.prior.scale"=current.disp.prior.scale.1)
  inits2 <- list("means0"=current.means0.2, "means1"=current.means1.2, 
                 "means2"=current.means2.2, "disps0"=current.disps0.2, 
                 "disps1"=current.disps1.2, "disps2"=current.disps2.2, 
                 "mean.prior.location"=current.mean.prior.location.2, 
                 "mean.prior.scale"=current.mean.prior.scale.2, 
                 "disp.prior.location"=current.disp.prior.location.2, 
                 "disp.prior.scale"=current.disp.prior.scale.2)
  inits3 <- list("means0"=current.means0.3, "means1"=current.means1.3, 
                 "means2"=current.means2.3, "disps0"=current.disps0.3, 
                 "disps1"=current.disps1.3, "disps2"=current.disps2.3, 
                 "mean.prior.location"=current.mean.prior.location.3, 
                 "mean.prior.scale"=current.mean.prior.scale.3, 
                 "disp.prior.location"=current.disp.prior.location.3, 
                 "disp.prior.scale"=current.disp.prior.scale.3)
  mean.proposal.scales0 <- current.mean.proposal.scales0
  mean.proposal.scales1 <- current.mean.proposal.scales1
  mean.proposal.scales2 <- current.mean.proposal.scales2
  disp.proposal.scales0 <- current.disp.proposal.scales0
  disp.proposal.scales1 <- current.disp.proposal.scales1
  disp.proposal.scales2 <- current.disp.proposal.scales2
  mean.prior.scale.proposal.sd <- current.mean.prior.scale.proposal.sd
  disp.prior.scale.proposal.sd <- current.disp.prior.scale.proposal.sd
  
  # Run MCMC ####
  message("running mcmc...")
  final.chains <- ln_hmm_3_chains(
    counts=counts, 
    groups=groups, 
    chain.length=chain.length, 
    thin=1, 
    inits1=inits1, 
    inits2=inits2, 
    inits3=inits3, 
    mean.proposal.scales0=mean.proposal.scales0, 
    mean.proposal.scales1=mean.proposal.scales1, 
    mean.proposal.scales2=mean.proposal.scales2, 
    disp.proposal.scales0=disp.proposal.scales0, 
    disp.proposal.scales1=disp.proposal.scales1, 
    disp.proposal.scales2=disp.proposal.scales2, 
    mean.prior.scale.proposal.sd=mean.prior.scale.proposal.sd, 
    disp.prior.scale.proposal.sd=disp.prior.scale.proposal.sd
  )
  rm(inits1, inits2, inits3)
  
  means0 <- mcmc.list(mcmc(final.chains$posterior.means0.1), 
                      mcmc(final.chains$posterior.means0.2), 
                      mcmc(final.chains$posterior.means0.3))
  means1 <- mcmc.list(mcmc(final.chains$posterior.means1.1), 
                      mcmc(final.chains$posterior.means1.2), 
                      mcmc(final.chains$posterior.means1.3))
  means2 <- mcmc.list(mcmc(final.chains$posterior.means2.1), 
                      mcmc(final.chains$posterior.means2.2), 
                      mcmc(final.chains$posterior.means2.3))
  disps0 <- mcmc.list(mcmc(final.chains$posterior.disps0.1), 
                      mcmc(final.chains$posterior.disps0.2), 
                      mcmc(final.chains$posterior.disps0.3))
  disps1 <- mcmc.list(mcmc(final.chains$posterior.disps1.1), 
                      mcmc(final.chains$posterior.disps1.2), 
                      mcmc(final.chains$posterior.disps1.3))
  disps2 <- mcmc.list(mcmc(final.chains$posterior.disps2.1), 
                      mcmc(final.chains$posterior.disps2.2), 
                      mcmc(final.chains$posterior.disps2.3))
  mean.prior.location <- mcmc.list(mcmc(final.chains$posterior.mean.prior.location.1), 
                                   mcmc(final.chains$posterior.mean.prior.location.2), 
                                   mcmc(final.chains$posterior.mean.prior.location.3))
  mean.prior.scale <- mcmc.list(mcmc(final.chains$posterior.mean.prior.scale.1), 
                                mcmc(final.chains$posterior.mean.prior.scale.2), 
                                mcmc(final.chains$posterior.mean.prior.scale.3))
  disp.prior.location <- mcmc.list(mcmc(final.chains$posterior.disp.prior.location.1), 
                                   mcmc(final.chains$posterior.disp.prior.location.2), 
                                   mcmc(final.chains$posterior.disp.prior.location.3))
  disp.prior.scale <- mcmc.list(mcmc(final.chains$posterior.disp.prior.scale.1), 
                                mcmc(final.chains$posterior.disp.prior.scale.2), 
                                mcmc(final.chains$posterior.disp.prior.scale.3))
  indicators <- mcmc.list(mcmc(final.chains$posterior.indicators.1), 
                          mcmc(final.chains$posterior.indicators.2), 
                          mcmc(final.chains$posterior.indicators.3))
  proportion <- mcmc.list(mcmc(final.chains$posterior.proportion.1), 
                          mcmc(final.chains$posterior.proportion.2), 
                          mcmc(final.chains$posterior.proportion.3))
  rm(final.chains)
  
  return(list("chain.length"=chain.length, 
              "adaptive.runs"=adapt.runs, 
              "mean.proposal.scales0"=mean.proposal.scales0, 
              "mean.proposal.scales1"=mean.proposal.scales1, 
              "mean.proposal.scales2"=mean.proposal.scales2, 
              "disp.proposal.scales0"=disp.proposal.scales0, 
              "disp.proposal.scales1"=disp.proposal.scales1, 
              "disp.proposal.scales2"=disp.proposal.scales2, 
              "mean.prior.scale.proposal.sd"=mean.prior.scale.proposal.sd, 
              "disp.prior.scale.proposal.sd"=disp.prior.scale.proposal.sd, 
              "means0"=means0, 
              "means1"=means1, 
              "means2"=means2, 
              "disps0"=disps0, 
              "disps1"=disps1, 
              "disps2"=disps2, 
              "mean.prior.location"=mean.prior.location, 
              "disp.prior.location"=disp.prior.location, 
              "mean.prior.scale"=mean.prior.scale, 
              "disp.prior.scale"=disp.prior.scale, 
              "indicators"=indicators, 
              "proportion"=proportion))
  
}
