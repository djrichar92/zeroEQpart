#' @title Calculate Confidence Interval of Equal Zero Order and (Semi) Partial
#'   Correlations

#' @description The \code{pzconf} function calculates confidence intervals for a
#'   zero order correlation minus a (semi) partial correlation (\code{p.xy -
#'   p.xy.z}). It is intended to be used with the \code{\link{pzcor}} function.
#'   
#' @param pzcor_obj pzcor object.
#'   
#' @param level numerical. Confidence level used to calculate the confidence 
#'   interval. This may be a vector so multiple intervals can be determined.
#'   
#' @details The \code{pzconf} function calculates confidence intervals based on 
#'   the bootstrap estimates determined from the \code{\link{pzcor}} function. 
#'   See \code{?pzcor} for details.
#'   
#' @return The confidence interval(s) is(are) displayed in a dataframe with four
#'   columns: Level, Lower, Upper, and Warnings. Level refers to the confidence 
#'   level of the interval. Lower and Upper are the respective lower and upper 
#'   bounds of the interval. Warnings may say "Max Level Passed" to show that 
#'   the specified confidence level is too high for the bootstrap estimate and 
#'   should be neglected. The last row (named "Max") is the largest confidence 
#'   level that can be determined by the given bootstrap estimate.
#'   
#' @seealso \code{\link{pzcor}}
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
#' test <- pzcor(data[,1], data[,2], data[,c(3,4)], k = 2000)
#' hist(test$distribution)
#' pzconf(test, c(0.9, 0.95, 0.99))
#'   
#' @export
pzconf = function(pzcor_obj, level = 0.9){
  # As per equation 14.10 on page 185 of "Introduction to the Bootstrap",
  # Efron and Tibshirani
  if(class(pzcor_obj) != 'pzcor'){
    first_arg = match.call()[2]
    stop('Must be pzcor object: ', first_arg)
  }
  z_0 <- pzcor_obj$bias
  acceleration <- pzcor_obj$acceleration
  distribution <- sort(pzcor_obj$distribution)
  k_eff <- pzcor_obj$k_eff
  test <- pzcor_obj$test
  
  max_level <- get_max_level(k_eff, acceleration, z_0, test)
  
  if(test == 'eq'){
    z_alpha <- stats::qnorm((1 - level)/2)
    result <- data.frame(Level = level, Lower = z_alpha, Upper = -z_alpha)
    max_result <- data.frame(Level = max_level, 
                             Lower = distribution[1], 
                             Upper = distribution[k_eff])
  } else {
    z_alpha <- stats::qnorm((1 - level))
    if(test == 'gt'){
      result <- data.frame(Level = level, Lower = z_alpha)
      max_result <- data.frame(Level = max_level, 
                               Lower = distribution[1])
    } else if(test == 'lt'){
      result <- data.frame(Level = level, Upper = -z_alpha)
      max_result <- data.frame(Level = max_level, 
                               Upper = distribution[k_eff])
    }
  }
  
  get_bound <- function(z_alpha){
    # As per equation 14.10 on page 185 of "Introduction to the Bootstrap",
    # Efron and Tibshirani
    BCa_z <- z_0 + (z_0 + z_alpha)/(1 - acceleration*(z_0 + z_alpha))
    cum_prop <- stats::pnorm(BCa_z)
    bound_ind <- cum_prop * k_eff
    bound_ind[which(bound_ind < 1)] <- 1
    lower_ind <- which(z_alpha < 0)
    upper_ind <- which(z_alpha >= 0)
    bound_ind[lower_ind] <- floor(bound_ind[lower_ind])
    bound_ind[upper_ind] <- ceiling(bound_ind[upper_ind])
    bound <- distribution[bound_ind]
    return(bound)
  }
  
  result[,-1] <- lapply(result[,-1, drop = FALSE], get_bound)
  
  result <- rbind(result, max_result)
  
  row.names(result)[nrow(result)] <- 'Max'
  
  result$Warnings <- ''
  
  result$Warnings[which(result$Level > max_level)] <- '  Max Level Passed'
  
  return(result)
}

get_max_level = function(k_eff, acceleration, bias, test){
  # As per equation 15.34 on page 216 of "Introduction to the Bootstrap",
  # Efron and Tibshirani
  min <- get_unadjusted_alpha(1/k_eff, acceleration, bias)
  max <- get_unadjusted_alpha(1 - 1/k_eff, acceleration, bias)
  if(test == 'eq'){
    result <- max - min
  } else if(test == 'gt'){
    result <- 1 - min
  } else if(test == 'lt'){
    result <- max
  }
  return(result)
}

get_unadjusted_alpha = function(cum_prop, acceleration, bias){
  w_0 <- stats::qnorm(cum_prop)
  z_0 <- bias
  BCa_z <- (w_0 - z_0)/(1 + acceleration*(w_0 - z_0)) - z_0
  alpha <- stats::pnorm(BCa_z)
  return(alpha)
  
}
