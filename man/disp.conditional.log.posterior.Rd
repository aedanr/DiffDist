% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conditional_posterior_functions_lognormal_hmm.R
\name{disp.conditional.log.posterior}
\alias{disp.conditional.log.posterior}
\title{Per-gene conditional log posterior dispersion function. Called by MCMC functions.}
\usage{
disp.conditional.log.posterior(
  genes,
  counts,
  n,
  sample.means,
  means,
  disps,
  prior.location,
  prior.scale
)
}
\arguments{
\item{genes}{Number of genes}

\item{counts}{Matrix of counts with genes in rows and samples in columns}

\item{n}{Number of samples}

\item{sample.means}{Vector of per-gene sample means}

\item{means}{Vector of current per-gene mean estimates}

\item{disps}{Vector of current per-gene dispersion estimates}

\item{prior.location}{Current estimate of prior location parameter}

\item{prior.scale}{Current estimate of prior scale parameter}
}
\value{
Vector of conditional posterior dispersion densities on log scale
}
\description{
Per-gene conditional log posterior dispersion function. Called by MCMC functions.
}
