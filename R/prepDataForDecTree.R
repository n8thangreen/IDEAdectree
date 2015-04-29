#' prepDataForDecTree
#'
#' \code{prepDataForDecTree} returns joined tables of inputs and predicted outcomes
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

prepDataForDecTree <- function(data, data.names){
  # relabel, remove and group variables
  #

  # removals

  othertest.names <- "IGRA"

  ## identify patient records with desired test outcomes (only)
  allTests <- apply(data[,c("TBcult","Smear")], 1, function(x)
    !any(x=="Not taken" | x=="INDETERMINATE" | x=="" | is.na(x)))

  # allTests <- allTests & apply(data[,c("TSTres","CXR")], 1, function(x) !any(x=="INDETERMINATE" | x==""))
#   allTests <- allTests & apply(data[,othertest.names, drop=FALSE], 1, function(x)
#       all(x=="Not taken" | x=="" | is.na(x))) # for when we don't want IGRA tested patients included

  ## IGRA or TST
  ## fudge!
  ## assuming that decisions from TST alone are the same as for IGRA or TST together then simply change column name
  ## remember to change column heading in ouput table...
  names(data)[names(data)=="TSTres"] <- "TSTres.old"
  names(data)[names(data)=="IGRAorTST"] <- "TSTres"

  ## replace 'INDETERMINATE' with 'Not taken', since assumed same because uninformative for diagnosis
  data[,c("TSTres","CXR")][t(apply(data[,c("TSTres","CXR")], 1, function(x) x%in%c("INDETERMINATE","")))] <- "Not taken"

  data <- data.frame(data[, data.names], Risk_factors=data$riskfacScore>0.6) ##TODO## check cut-off value
  data <- data[allTests,]

  data$treatResponse[is.na(data$treatResponse)] <- FALSE

  data
}
