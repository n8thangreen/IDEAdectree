
#' Make rule-out table for multiple follow-up times
#'
#' Iterate over a range of follow-up times and row-bind results (potentially a high number of rows).
#'
#' @param FNtime.seq
#' @param Ctest
#'
#' @return out

make.ruleoutTable.FNtime <- function(FNtime.seq=c(7,seq(0,60,by=1)), Ctest){
  out <- NULL
  for (i in FNtime.seq){
    out <- rbind(out,
                 cbind(make.ruleoutTable.pre(FNtime=i, Ctest=Ctest,
                                             thresh = seq(from=1, to=0.8, by=-0.005))$combinedDosanjh, FNtime=i))
  }
  out
}


#' Make rule-out table of multiple clinically judged highrisk cut-offs
#'
#' Iterate over a range of cut-off values [0,1] and row-bind results (potentially a high number of rows).
#'
#' @param highrisk.seq
#' @param Ctest
#'
#' @return out

make.ruleoutTable.highrisk <- function(highrisk.seq=seq(0.1,0.99,by=0.01), Ctest){
  out <- NULL
  for (i in highrisk.seq){
    out <- rbind(out,
                 cbind(make.ruleoutTable.pre(prop_highrisk=i, Ctest=Ctest,
                                             thresh = seq(from=1, to=0.8, by=-0.005))$combinedDosanjh, prop_highrisk=i))
  }
  out
}


#' Make rule-out table for multiple cohort prevalences
#'
#' Iterate over a range of prevalences 0-100% and row-bind results (potentially a high number of rows).
#'
#' @param cat4percent.seq
#' @param Ctest
#'
#' @return out

make.ruleoutTable.cat4percent <- function(cat4percent.seq = seq(1,100,by=1), Ctest){
  out <- NULL
  for (i in cat4percent.seq){
    out <- rbind(out,
                 cbind(make.ruleoutTable.pre(cat4percent=i, Ctest=Ctest,
                                             thresh = seq(from=1, to=0.8, by=-0.005))$combinedDosanjh, cat4percent=i))
  }
  out
}


#' Make rule-out table of bootstrapped samples of patients
#'
#' Iterate over bootstrap samples and row-bind results (potentially a high number of rows).
#'
#' @param nboot
#' @param Ctest
#'
#' @return out

make.ruleoutTable.bootsample <- function(nboot = 2, Ctest){
  out <- NULL
  for (i in 1:nboot){

    bootrows <- sample(1:nrow(data), nrow(data), replace=TRUE)
    bootdata <- data[bootrows,]
    out <- rbind(out,
                 cbind(make.ruleoutTable.pre(data=bootdata, #need to change original function!
                                             Ctest=Ctest, thresh=seq(from=1, to=0.5, by=-0.05))$combinedDosanjh, bootsample=i))
  }
  out
}

