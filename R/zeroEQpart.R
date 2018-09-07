#' @title Test for Equal Zero Order and (Semi) Partial Correlations

#' @description Calculate the statistical significance of a zero order
#'   correlation being equal to a partial or semi-partial correlation using
#'   the bias-corrected and accelerated bootstrap method.
#'   
#' @section pzcor: The \code{pzcor} function performs one of the specified 
#'   tests: \itemize{ \item{\code{p.xy - p.xy.z = 0} (default)} \item{\code{p.xy
#'   - p.xy.z > 0}} \item{ \code{p.xy - p.xy.z < 0}} } See \code{?pzcor} for 
#'   details.
#'   
#' @section pzconf: The \code{pzconf} function computes a confidence interval 
#'   for the statistic: \code{p.xy - p.xy.z}. To be used with \code{pzcor}. See 
#'   \code{?pzconf} for details.
#'   
#' @docType package
#' @name zeroEQpart
#'   
NULL
