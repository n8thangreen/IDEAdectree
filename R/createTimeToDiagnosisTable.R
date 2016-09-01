
#' Get All Data Test Combinations
#'
#' \code{get.allDataTestCombinations}
#'
#' @param data
#' @return data
#'

get.allDataTestCombinations <- function(data){

  test.names <- c("IGRA", "Smear", "TBcult",# "QFN", "TSPOT", "CSF", "BAL",
                  "HistBiop", "NeedleAsp", "PCR", "imaging")# "PET", "CXR", "CT", "MRI")
  pathdata <- data[,test.names]
  names(pathdata) <- paste(names(pathdata), "bin", sep="_")
  pathdata <- as.matrix(pathdata)
  pathdata[pathdata!="Not taken"] <- 1
  pathdata[pathdata=="Not taken"] <- 0
  data <- data.frame(data, pathdata)
  pathdata <- as.data.frame(pathdata[!duplicated(pathdata),])
  pathdata <- pathdata[do.call(order, pathdata),]
  data.frame(pathID=1:nrow(pathdata), pathdata)
}


#' get.statStartToDiagByCovariate
#'
#' \code{get.statStartToDiagByCovariate}
#'
#' @param data
#' @return data
#'

get.statStartToDiagByCovariate <- function(field, value, stat){
  aggregate(start.to.diag ~ pathID, data=data.joined[data.joined[,field]%in%value,], stat, na.rm=T)
}


#' get.meansdStartToDiagByCovariate
#'
#' \code{get.meansdStartToDiagByCovariate}
#'
#' @param data
#' @return data
#'

get.meansdStartToDiagByCovariate <- function(field, value, cutoff=8){

  xbar <- get.statStartToDiagByCovariate(field, value, mean)
  stdev <- get.statStartToDiagByCovariate(field, value, sd)
  x <- merge(xbar, stdev, by="pathID", all=TRUE)
  names(x)[2] <- paste(value, "mean", sep="_")
  names(x)[3] <- paste(value, "sd", sep="_")

  quants <- as.matrix(get.statStartToDiagByCovariate(field, value, quantile))
  colnames(quants) <- c("pathID", paste(c("0", "25", "50", "75", "100"), value, sep="_"))
  x <- merge(x, quants, by="pathID", all=TRUE)
  names(x) <- sub("start.to.diag", value, names(x))

  n <- aggregate(start.to.diag ~ pathID, data=data.joined[data.joined[,field]==value,], length)
  x <- merge(x, n, by="pathID", all=TRUE)
  names(x)[length(names(x))] <- paste("N", value, sep="_")

  x <- x[x$N>cutoff,]

  return(x)
}
