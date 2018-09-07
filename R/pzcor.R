#' @title Calculate P-Value of Equal Zero Order and (Semi) Partial Correlations

#'@description Calculate the p-value of a zero order correlation being equal to
#'  a partial or semi-partial correlation using bootstrap.
#'@param x a numeric vector.
#'
#'@param y a numeric vector.
#'
#'@param z a numeric vector.
#'
#'@param semi logical. If \code{TRUE}, then the semi-partial correlation between
#'  \code{x} and \code{y} given \code{z} is used. If \code{FALSE} (default),
#'  then the partial correlation between \code{x} given \code{z} and \code{y}
#'  given \code{z} is used.
#'
#'@param k the number of bootstrap samples taken (default is 1000).
#'
#'@param method a character string indicating which partial correlation
#'  coefficient is to be computed. One of "pearson" (default), "kendall", or
#'  "spearman" can be abbreviated.
#'
#'@param test character string denoting the type of test to be administered. Can
#'  be one of the three: \itemize{ \item{\code{'eq'} tests  \code{p.xy - p.xy.z =
#'  0} (default)} \item{\code{'gt'} tests  \code{p.xy - p.xy.z > 0}}
#'  \item{\code{'lt'} tests  \code{p.xy - p.xy.z < 0}} }
#'
#'@details Uses the bias-corrected and accelerated (BCa) bootstrap method to
#'  estimate the distribution of the difference \eqn{\theta = \rho_xy -
#'  \rho_xy|z} where \eqn{\rho_xy} is the zero order correlation between
#'  variables \eqn{x} and \eqn{y} and \eqn{\rho_xy|z} is the (semi) partial
#'  correlation between the respective variables while accounting for those in
#'  \eqn{z}. If \eqn{|\theta| > 0}, then the p-value will not be calculated.
#'
#'@return \item{acceleration}{the acceleration used for the BCa method.}
#'
#'  \item{alpha}{the proportion of the bootstrapped distribution below zero.}
#'
#'  \item{bias}{the bias used for the BCa method.}
#'
#'  \item{call}{shows the function call.}
#'
#'  \item{distribution}{the estimated distribution of the difference as
#'  determined through bootstrapping.}
#'
#'  \item{k_eff}{the number of successful bootstrap samples. Less than or equal
#'  to \code{k}.}
#'
#'  \item{method}{the method of correlation used.}
#'
#'  \item{p.value}{significance level of the test.}
#'
#'  \item{p.xy}{Zero order correlation between \code{x} and \code{y}.}
#'
#'  \item{p.xy.z}{(semi) partial correlation between \code{x} and \code{y} while
#'  accounting for \code{z}.}
#'
#'  \item{semi}{logical. If \code{TRUE}, \code{p.xy.z} is the semi-partial
#'  correlation. Otherwise \code{p.xy.z} is the partial correlation.}
#'
#'  \item{difference}{calculated from the data. Same as \code{p.xy - p.xy.z}.}
#'
#'  \item{test}{shows the type of test performed.}
#'
#'@seealso \code{\link{pzconf}}
#'
#' @examples
#' require(graphics)
#' require(MASS)
#' # data
#' set.seed(1111)
#' mu <- rep(0,4)
#' Sigma <- matrix(.2, nrow=4, ncol=4) + diag(4)*.8
#' data <- mvrnorm(n=100, mu=mu, Sigma=Sigma)
#'
#' # p.(1,2) = p.(1,2)|(3,4) test
#' test <- pzcor(data[,1], data[,2], data[,c(3,4)], k = 2000, semi = FALSE, test = 'eq')
#' hist(test$distribution)
#' test

#' @export
pzcor <- function(x, y, z, semi = FALSE, k = 1000, method = "pearson",
                  test = 'eq'){
  # generic function for S3 class: pzcor
  UseMethod("pzcor")
}

#' @export
pzcor.default <- function(x, y, z, semi = FALSE, k = 1000, method = "pearson",
                          test = 'eq'){
  # Default method for pzcor
  # Calls various helper functions and returns the results in a list

  method <- get_cor_method(method)

  dist_result <- get_dist(x,y,z,semi,k,method)

  alpha <- get_alpha(dist_result$distribution)

  bias <- get_bias(dist_result$distribution, dist_result$difference)

  acceleration <- get_acceleration(dist_result$jack_distribution)

  k_eff <- dist_result$k_eff

  p.value <- get_p.value(alpha, acceleration, bias, test, k_eff)

  result <- list(
                acceleration = acceleration,

                alpha = alpha,

                bias = bias,

                call = match.call(),

                distribution = dist_result$distribution,

                k_eff = k_eff,

                method = method,

                p.value = p.value,

                p.xy = dist_result$p.xy,

                p.xy.z = dist_result$p.xy.z,

                semi = semi,

                difference = dist_result$difference,

                test = test
                )

  class(result) <- "pzcor"

  return(result)
}

#' @export
print.pzcor <- function(x, ...){
  # print method for pzcor objects
  print(summary(x))

}

#' @export
summary.pzcor <- function(object, ...){
  if(object$test == 'eq'){
    hyp_op <- '='
  } else if(object$test == 'gt'){
    hyp_op <- '>'
  } else if(object$test == 'lt'){
    hyp_op <- '<'
  }

  hypothesis <- paste('p.xy - p.xy.z ', hyp_op, ' 0')
  test_summary <- list(
                            'Hypothesis               ' = hypothesis,

                            Semi_Partial = object$semi,

                            Correlation = object$method,

                            p.xy = object$p.xy,

                            p.xy.z = object$p.xy.z,

                            difference = object$difference,

                            p.value = object$p.value,

                            Bootstrap_Size = object$k_eff
  )

  relation <- 'equal to'
  if(object$alpha %in% c(0,1)){
    if(object$p.value < 0.5){
      relation <- 'less than'
    } else {
      relation <- 'greater than'
    }
  }

  pval_name <- paste('p.value ', relation)

  names(test_summary)[which(names(test_summary)=='p.value')] <- pval_name

  distr_summary <- summary(object$distribution)

  result <- list(test_summary = test_summary,
                 distribution_summary = distr_summary,
                 call = object$call)

  class(result) <- 'pzcor_summary'

  return(result)
}

#' @export
print.pzcor_summary <- function(x, ...){

  # put test_summary in proper format
  test_summary <- x$test_summary

  b_size <- test_summary$Bootstrap_Size

  is.num <- sapply(rbind(test_summary), is.numeric)

  form_nums <- formatC(unlist(test_summary[is.num]), digits = 5, format = "f")

  test_summary[is.num] <- form_nums

  test_summary$Bootstrap_Size <- b_size

  test_summary_formatted <- data.frame(cbind(test_summary))

  colnames(test_summary_formatted) <- ''


  # print summary
  print(x$call)
  stars <- rep('*', 8, sep='')
  cat('\n', stars, 'Test Summary', stars)

  print(test_summary_formatted)

  stars <- rep('*', 9, sep='')
  cat('\n\n', stars, 'Distribution Summary', stars, '*', '\n')

  print(x$distribution_summary)
}

get_dist <- function(x, y, z, semi = FALSE, k = 1000, method = "pearson"){
  # Estimates distribution function of theta = p_{xy} - p_{xy|z} (or p_{xy} -
  # p_{x, y|z} if semi is TRUE) using bootstrap and jackknife (for
  # acceleration). Returns a list including both distributions, p_{xy}, p_{xy|z}
  # (or p_{x, y|z}), and theta.
  # Used as a helper function for pzcor

  x <- data.matrix(x)
  y <- data.matrix(y)
  z <- data.matrix(z)

  zero_ord_cor <- function(ind){
    p.xy <- stats::cor(x[ind], y[ind], method = method)
    return(p.xy)
  }

  part_cor <- function(ind){
    p.xy.z <- ppcor::pcor.test(x[ind], y[ind], z[ind,], method = method)
    return(p.xy.z$estimate)
  }

  semi_part_cor <- function(ind){
    p.xy.z <- ppcor::spcor.test(x[ind], y[ind], z[ind,], method = method)
    return(p.xy.z$estimate)
  }

  n <- length(x)

  boot_ind  <- data.matrix(replicate(k, sample(1:n, replace = TRUE)))
  jack_ind  <- data.matrix(utils::combn(x = 1:n, m = n - 1))

  difference_p.xy <- zero_ord_cor(1:n)
  p.xy <- apply(X = boot_ind, FUN = zero_ord_cor, MARGIN = 2)

  jack_p.xy <- apply(X = jack_ind, FUN = zero_ord_cor, MARGIN = 2)

  if (semi == TRUE){
    difference_p.xy.z <- semi_part_cor(1:n)
    p.xy.z <- apply(X = boot_ind, FUN = semi_part_cor, MARGIN = 2)
    jack_p.xy.z <- apply(X = jack_ind, FUN = semi_part_cor, MARGIN = 2)
  } else {
    difference_p.xy.z <- part_cor(1:n)
    p.xy.z <- apply(X = boot_ind, FUN = part_cor, MARGIN = 2)
    jack_p.xy.z <- apply(X = jack_ind, FUN = part_cor, MARGIN = 2)
  }

  boot_data = stats::na.omit(data.frame(p.xy, p.xy.z))
  jack_data = stats::na.omit(data.frame(jack_p.xy, jack_p.xy.z))


  result <- list()

  result$distribution <- boot_data[,1] - boot_data[,2]

  result$jack_distribution <- jack_data[,1] - jack_data[,2]

  result$k_eff <- nrow(boot_data)

  result$p.xy <- difference_p.xy

  result$p.xy.z <- difference_p.xy.z

  result$difference <- difference_p.xy - difference_p.xy.z

  return(result)
}

p_less_x <- function(dist, x){
  # Returns the proportion of values in dist that are less than x
  n <- length(dist)

  less_than_x <- length(dist[dist < x])
  p_x <- less_than_x/n

  return(p_x)
}

get_alpha <- function(dist){
  # As per equation 15.31 on page 215 of "Introduction to the Bootstrap",
  # Efron and Tibshirani
  alpha = p_less_x(dist, 0)
  return(alpha)
}

get_acceleration <- function(dist){
  # Receives jackknife distribution and determines acceleration for p.value
  # adjustment.
  # As per equation 14.15 on page 186 of "Introduction to the Bootstrap",
  # Efron and Tibshirani
  x_dist <- mean(dist) - dist
  num <- sum(x_dist^3)
  denom <- 6*sum(x_dist^2)^(3/2)
  return(num/denom)

}

get_bias <- function(distribution, difference){
  # As per equation 14.14 on page 186 of "Introduction to the Bootstrap",
  # Efron and Tibshirani
  bias_proportion <- p_less_x(distribution, difference)
  bias_z <- stats::qnorm(bias_proportion)
  return(bias_z)
}

get_p.value <- function(alpha, acceleration, bias, test, k_eff){
  # As per equation 15.34 on page 216 of "Introduction to the Bootstrap",
  # Efron and Tibshirani
  if(alpha == 0){
    alpha <- 1/k_eff
  } else if(alpha == 1){
    alpha <- 1 - 1/k_eff
  }

  w_0 <- stats::qnorm(alpha)
  z_0 <- bias

  BCa_z <- (w_0 - z_0)/(1 + acceleration*(w_0 - z_0)) - z_0
  p.value <- 2*stats::pnorm(-1*abs(BCa_z))

  if(test == 'gt' & BCa_z > 0 | test == 'lt' & BCa_z < 0){
      p.value <- p.value/2

  } else if(test != 'eq'){
    p.value <- 1 - p.value/2
  }

  return(p.value)
}

get_cor_method <- function(abr_method){
  cors <- c('pearson', 'kendall', 'spearman')
  method <- match.arg(abr_method, cors)
  return(method)
}

