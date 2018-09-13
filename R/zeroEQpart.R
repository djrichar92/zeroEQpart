#' @title Zero Order vs (Semi) Partial Correlation Test and CI
#'
#' @description Calculate the statistical significance of a zero order
#'   correlation being equal to a partial or semi-partial correlation using
#'   the bias-corrected and accelerated bootstrap method from "An Introduction
#'   to the Bootstrap" Efron (1983) <0-412-04231-2>. Confidence intervals for
#'   the parameter (zero order minus partial) can also be determined.
#'
#' @section pzcor: The \code{pzcor} function tests one of the following null
#'   hypotheses: \itemize{
#'                   \item{\eqn{\rho.xy - \rho.xy.z = 0} (default)}
#'                   \item{\eqn{\rho.xy - \rho.xy.z \ge 0}}
#'                   \item{ \eqn{\rho.xy - \rho.xy.z \le 0}}
#'                   }
#'   See \code{\link{pzcor}} for details.
#'
#' @section pzconf: The \code{pzconf} function computes confidence intervals
#'   for the parameter: \eqn{\rho.xy - \rho.xy.z}. To be used with
#'   \code{pzcor}. See \code{\link{pzconf}} for details.
#'
"_PACKAGE"
