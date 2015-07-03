#' Rule-based time to diagnosis estimation
#'
#' \code{estimateTimeToDiagnosis} estimates time to diagnosis for each patient using simple pre-defined rules.
#' This time is not directly available in the data but can be got at in most cases indirectly.
#'
#' @param data Individual patient records
#' @return data with the start.to.diag column appended


estimateTimeToDiagnosis <- function(data){
  ## call: data <- estimateTimeToDiagnosis(data)

  ## the rules for estimationg time to diagnosis are:
  ## 1) If culture +ve then take report time
  ## 2) What is means of diagnosis? Use corresponding time.
  ## 3) min{DateDiagCon, lasttestDate}


  is.TestUsedAsMeansOfDiag <- function(data, i, testKeyword){

    datai <- data[i, c("MeansFinDiag1", "MeansFinDiagother.multiple", "Diagnostic_mean")]
    any(sapply(datai, function(x) grepl(testKeyword, x, ignore.case=TRUE)))
  }


  getDiagnosisTimeViaMeans <- function(data, i){

    diagtime <- NA
    diag.list <- list("Smear"="start.to.Smear", "Imaging"="start.to.Imaging", "Histology"="start.to.HistBiop",
                      "IGRA"="start.to.IGRA", "Culture"="start.to.TBcult", "Response to Treatment"="testDrug_diff_plus63days",
                      "PET"="start.to.PET", "BAL"="start.to.BAL", "PCR"="start.to.PCR", "TST"="start.to.TST")
                      # "EBUS"=NA, "Clinical features"=NA, "Empiric"=NA)  #as date of presentation?

    for (testKeyword in names(diag.list)){

      if(is.TestUsedAsMeansOfDiag(data, i, testKeyword)) diagtime <- c(diagtime, data[i, diag.list[[testKeyword]]])
    }
    max(diagtime, na.rm=TRUE) #min?
  }


  for (i in 1:nrow(data)){

    if (data$TBcult[i]=="POSITIVE")
      data$start.to.diag[i]  <- data$testCult_diff[i]
    else{
      data$start.to.diag[i] <- getDiagnosisTimeViaMeans(data, i)}
  }

  return(data)
}



