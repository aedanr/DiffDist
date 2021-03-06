% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conditional_posterior_functions_lognormal_hmm.R
\name{mean.prior.scale.log.posterior}
\alias{mean.prior.scale.log.posterior}
\title{Conditional log posterior function for scale parameter of mean prior. Called by
MCMC functions.}
\usage{
\method{mean}{prior.scale.log.posterior}(prior.scale, prior.location, g, means0, means1, means2, z)
}
\arguments{
\item{prior.scale}{Current estimate of prior scale parameter}

\item{prior.location}{Current estimate of prior location parameter}

\item{g}{Number of genes}

\item{means0}{Vector of current overall per-gene mean estimates}

\item{means1}{Vector of current group 1 per-gene mean estimates}

\item{means2}{Vector of current group 2 per-gene mean estimates}

\item{z}{Vector of current per-gene mixture component indicator variables}
}
\value{
Conditional posterior density of scale parameter of mean prior, on log scale
}
\description{
Conditional log posterior function for scale parameter of mean prior. Called by
MCMC functions.
}
