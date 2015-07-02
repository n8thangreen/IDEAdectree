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


  for (i in 1:nrow(data)){
    eachlist <- vector("numeric", 0)

    if (data$TBcult[i]=="POSITIVE")
      eachlist <- c(testCult_diff=data$testCult_diff[i])
    else{

      diagcols.names <- c("MeansFinDiag1", "MeansFinDiagother.multiple", "Diagnostic_mean")
      diag.list <- list("Smear"=start.to.Smear, "Imaging"=start.to.Imaging, "Histology"=start.to.HistBiop, "Clinical features"=NA,
                        "IGRA"=start.to.IGRA, "Culture"=start.to.TBcult, "Response to Treatment"=TBDrugStart.min,
                        "EBUS"=NA, "PET"=start.to.PET, "BAL"=start.to.BAL, "PCR"=start.to.PCR, "TST"=start.to.TST, "Empiric"=NA)

      for (keyword in names(diag.list)){

        testTF <- any(sapply(data[i, diagcols.names], function(x) grepl(keyword, x, ignore.case=TRUE)))
        if(testTF)  diagtime <- diag.list[[keyword]]
      }
      data$start.to.diag[i] <- data[i,diagtime]
      eachlist <- c(eachlist, data[i,diagtime])

    }


    alllist[[i]] <- eachlist
  }

  # lapply(alllist, function(x) names(which(x==min(x, na.rm=T), arr.ind=TRUE)))

  combinedtimes <- unlist(lapply(alllist, min, na.rm=T))
  combinedtimes[is.infinite(combinedtimes)] <- NA

  return(round(combinedtimes))

}



