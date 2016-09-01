#' Rule-based Time to Diagnosis Estimation
#'
#' \code{estimateTimeToDiagnosis} estimates time to diagnosis for each patient using simple pre-defined rules.
#' This time is not directly available in the data but can be got at in most cases indirectly.
#'
#' The rules for estimating time to diagnosis are in priority order:
#' \enumerate{
#'  \item {If culture +ve then \code{min}(report time, start treatment) [previously, take report time only.]}
#'  \item {What is means of diagnosis? Use corresponding times. Latest date of those for the means of diagnosis used.}
#'  \item If diagnosis date is too late (130 days after first test) then use \code{DateDiagCon} if <130 days
#'  \item If diagnosis date still \code{NA} then use  \code{DateDiagCon}.
#'  \item {could use \code{min(DateDiagCon, lasttestDate)} isntead?}
#'  }
#'
#' @param data Individual patient records
#' @return data with the \code{start.to.diag} column appended
#' @export
estimateTimeToDiagnosis <- function(data){

  is.TestUsedAsMeansOfDiag <- function(data, i, testKeyword){

    datai <- data[i, c("MeansFinDiag1", "MeansFinDiagother.multiple", "Diagnostic_mean")]
    any(sapply(datai, function(x) grepl(testKeyword, x, ignore.case=TRUE)))
  }


  getDiagnosisTimeViaMeans <- function(data, i){

    diagtime  <- NA
    diag.list <- list("Smear" = "start.to.Smear",
                      "Imaging" = "start.to.Imaging",
                      "Histology" = "start.to.Histology",
                      "IGRA" = "start.to.IGRA",
                      "Culture" = "start.to.TBcultorig",
                      "Response to Treatment" = "start.to.FU",
                      "OTHER" = "start.to.other",
                      "PET" = "start.to.PET",
                      "BAL" = "start.to.BAL",
                      "PCR" = "start.to.PCR",
                      "TST" = "start.to.TST",#)
                      "Clinical features" = "start.to.clinicalfeatures") #as date of presentation?
    # "EBUS"=NA, "Empiric"=NA)

    for (testKeyword in names(diag.list)){

      if(is.TestUsedAsMeansOfDiag(data, i, testKeyword)) diagtime <- c(diagtime, data[i, diag.list[[testKeyword]]])
    }
    return(max(diagtime, na.rm=TRUE))
  }

  for (i in 1:nrow(data)){

    if (data$TBcult[i]=="POSITIVE"){
      data$start.to.diag[i] <- min(data$testDrug_diff[i], data$start.to.TBcultorig[i], na.rm = TRUE)
      data$start.to.diag[i] <- ifelse(is.na(data$start.to.diag[i]), 0, data$start.to.diag[i])

    }else{
      data$start.to.diag[i] <- getDiagnosisTimeViaMeans(data, i)}
  }

  data$start.to.diag[is.infinite(data$start.to.diag)] <- NA
  data$start.to.diag[data$start.to.diag<0] <- NA
  data$testDiagCon_diff[data$testDiagCon_diff<0] <- NA

  tooBigLimit <- 130  #days
  tooBig <- data$start.to.diag>tooBigLimit & data$testDiagCon_diff<tooBigLimit &
    !is.na(data$start.to.diag) & !is.na(data$testDiagCon_diff)
  data$start.to.diag[tooBig] <- data$testDiagCon_diff[tooBig]
  data$start.to.diag[is.na(data$start.to.diag)] <- data$testDiagCon_diff[is.na(data$start.to.diag)]

  invisible(data)
}

