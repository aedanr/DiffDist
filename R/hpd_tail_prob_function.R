#' Highest density interval function, used in \code{hpd.pval()}.
#'
#' Based on \code{HDInterval::hdi()} package, but with sorting step removed since
#' sorting is done in \code{hpd.pval()}.
#' @param x (Posterior) sample
#' @param credMass Desired density within interval
#' @return (lower, upper) boundaries of the narrowest interval containing
#' the desired density
#' @export
hdi <- function(x, credMass) {
  n <- length(x)
  exclude <- n - floor(n * credMass)
  low.poss <- x[1:exclude]
  upp.poss <- x[(n - exclude + 1):n]
  best <- which.min(upp.poss - low.poss)
  result <- c(low.poss[best], upp.poss[best])
  return(result)
}

#' Function to calculate posterior tail probability that \code{x > m} using highest
#' density intervals.
#' If \code{x} is a posterior sample of the log-transformed difference in a given
#' parameter between groups, \code{hpd.pval()} estimates the posterior probability that
#' the log fold change in the parameter between groups is greater than the threshold
#' \code{m}.
#' @param x Posterior sample
#' @param m Threshold probability. Default is 0.
#' @return Posterior probability that \code{x > m}
#' @export
hpd.pval <- function(x, m=0) {
  x <- sort.int(x, method='quick')

  y1 <- numeric(9)
  for (t in 1:9) {
    y1[t] <- sign(hdi(x, credMass=1 - t * 0.1)[1] - m) ==
      sign(hdi(x, credMass=1 - t * 0.1)[2] + m)
    if (y1[t] == 1) {break}
  }
  if (sum(y1) == 0) {x1 <- 9} else {x1 <- min(which(y1 == 1)) - 1}
  iter <- x1 * 0.1

  y2 <- numeric(9)
  for (t in 1:9) {
    y2[t] <- sign(hdi(x, credMass=1 - (iter + t * 0.01))[1] - m) ==
      sign(hdi(x, credMass=1 - (iter + t * 0.01))[2] + m)
    if (y2[t] == 1) {break}
  }
  if (sum(y2) == 0) {x2 <- 9} else {x2 <- min(which(y2 == 1)) - 1}
  iter <- iter + x2 * 0.01

  y3 <- numeric(9)
  for (t in 1:9) {
    y3[t] <- sign(hdi(x, credMass=1 - (iter + t * 1e-3))[1] - m) ==
      sign(hdi(x, credMass=1 - (iter + t * 1e-3))[2] + m)
    if (y3[t] == 1) {break}
  }
  if (sum(y3) == 0) {x3 <- 9} else {x3 <- min(which(y3 == 1)) - 1}
  iter <- iter + x3 * 1e-3

  y4 <- numeric(9)
  for (t in 1:9) {
    y4[t] <- sign(hdi(x, credMass=1 - (iter + t * 1e-4))[1] - m) ==
      sign(hdi(x, credMass=1 - (iter + t*1e-4))[2] + m)
    if (y4[t] == 1) {break}
  }
  if (length(x) <= 1e4) {
    if (sum(y4) == 0) {x4 <- 10} else {x4 <- min(which(y4 == 1))}
    iter <- iter + x4 * 1e-4
  }
  else {
    if (sum(y4) == 0) {x4 <- 9} else {x4 <- min(which(y4 == 1)) - 1}
    iter <- iter + x4 * 1e-4

    y5 <- numeric(9)
    for (t in 1:9) {
      y5[t] <- sign(hdi(x, credMass=1 - (iter + t * 1e-5))[1] - m) ==
        sign(hdi(x, credMass=1 - (iter + t * 1e-5))[2] + m)
      if (y5[t] == 1) {break}
    }
    if (sum(y5) == 0) {x5 <- 10} else {x5 <- min(which(y5 == 1))}
    iter <- iter + x5 * 1e-5
  }

  return(iter)
}
