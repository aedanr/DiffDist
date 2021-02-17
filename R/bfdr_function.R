#' Bayesian FDR function, based on ShrinkBayes::BFDR.
#' @param x Posterior probabilities
#' @return Corresponding Bayesian FDRs
#' @export
bfdr <- function(x) {
  ord <- sort(x, decreasing=T, index.return=T)
  bfdr <- cumsum(1-ord$x)/(1:length(x))
  bfdr[ord$ix] <- bfdr
  return(bfdr)
}
