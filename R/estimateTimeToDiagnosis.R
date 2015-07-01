#' Rule-based time to diagnosis estimation
#'
#' \code{estimateTimeToDiagnosis} estimates time to diagnosis for each patient using simple pre-defined rules
#'
#' @param data Individual patient records
#' @return timeToDiag

estimateTimeToDiagnosis <- function(data){
  ## call: estimateTimeToDiagnosis(data)

  ## the rules for estimationg time to diagnosis are:
  ## 1) If culture +ve then take report time
  ## 2) What is means of diagnosis? Use corresponding time.
  ## 3) min{DateDiagCon, lasttestDate}



  isValid <- function(x) (x>=0) & !is.na(x)

  alllist <- list()

  data$testCultorig_diff <- difftime(data$TBculttestDate.orig, data$testDate.min, units="days")
  testDrugEnd <- data$testDrug_diff + data$TBDrug_diff
  testDate.names <- c("TBculttestDate","QFNtestDate", "TSPOTtestDate", "TSTtestDate", "SmeartestDate", "BALtestDate",
                      "HistBioptestDate", "NeedleAsptestDate", "PCRtestDate", "CXRtestDate", "CTtestDate",
                      "MRItestDate", "Othertest1Date",  "Othertest2Date","PETtestDate")
  data$testDate.max <- apply(data[,testDate.names[-1]], 1, max, na.rm=T)
  data$testminmax_diff <- difftime(data$testDate.max, data$testDate.min, units="days")

  for (i in 1:nrow(data)){
    eachlist <- vector("numeric", 0)
    if (data$TBcult[i]=="POSITIVE")
      eachlist <- c(testCult_diff=data$testCult_diff[i])
    else{
      if (isValid(data$testCultorig_diff[i]) & (data$TBcult[i]=="NEGATIVE") & !data$TBcultCens[i])
        eachlist <- c(testCultorig_diff=data$testCultorig_diff[i])

      if (isValid(data$testDrug_diff[i]) & (data$Dosanjh[i]=="2"))
        eachlist <- c(eachlist, testDrug_diff=data$testDrug_diff[i] + 62)

      if (isValid(testDrugEnd[i]) & (data$DosanjhGrouped[i]=="4"))
        eachlist <- c(eachlist, testDrugEnd=data$testDrugEnd[i])

      if (isValid(data$testminmax_diff[i]) & (data$DosanjhGrouped[i]=="4"))
        eachlist <- c(eachlist, testminmax_diff=data$testminmax_diff[i])
    }
    alllist[[i]] <- eachlist
  }

  # lapply(alllist, function(x) names(which(x==min(x, na.rm=T), arr.ind=TRUE)))

  combinedtimes <- unlist(lapply(alllist, min, na.rm=T))
  combinedtimes[is.infinite(combinedtimes)] <- NA


  return(round(combinedtimes))
}



