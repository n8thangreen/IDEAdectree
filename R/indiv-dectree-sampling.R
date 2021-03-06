
#' Get standard deviation from Normal confidence interval
#'
#' @param n Sample size
#' @param x_bar Mean
#' @param upperCI Upper 95% CI
#' @param lowerCI Lower 95% CI
#'
#' @return sd

get.sd.from.normalCI <- function(n, x_bar=NA, upperCI=NA, lowerCI=NA){
  # http://stats.stackexchange.com/questions/30402/how-to-calculate-mean-and-standard-deviation-in-r-given-confidence-interval-and

  if(!is.na(lowerCI) & !is.na(x_bar)){
    sd <- sqrt(n) * (x_bar - lowerCI)/1.96
  }else if(!is.na(upperCI) & !is.na(x_bar)){
    sd <- sqrt(n) * (upperCI - x_bar)/1.96
  }else if (!is.na(lowerCI) & !is.na(upperCI)){
    sd <- sqrt(n) * (upperCI - lowerCI)/(2*1.96)
  }
  return(sd)
}


#' Method of Moments Beta Distribution Parameter Transformation
#'
#' Could alternatively use the beta-PERT \code{\link{rpert}} with maximum and minimum instead of variance.
#'
#' @param xbar Mean
#' @param vbar Variance
#'
#' @return a and b of Beta(a,b)
#' @seealso rpert

MoM.beta <- function(xbar, vbar){

  if(vbar==0){stop("zero variance not allowed")
  }else if(xbar*(1-xbar)<vbar){
    stop("mean or var inappropriate")
  }else{
    a <- xbar * (((xbar*(1-xbar))/vbar)-1)
    b <- (1-xbar) * (((xbar*(1-xbar))/vbar)-1)
  }
  list(a=a, b=b)
}


#' Method of Moments Gamma Distribution Parameter Transformation
#'
#' @param mean Mean
#' @param var Variance
#'
#' @return shape, scale
#' @seealso MoM.beta

MoM.gamma <- function(mean, var){

  stopifnot(var>=0)
  stopifnot(mean>=0)
  names(mean) <- NULL
  names(var)  <- NULL

  list(shape = mean^2/var,
       scale = var/mean)
}
