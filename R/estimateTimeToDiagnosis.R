
#' Rules-based Time to Diagnosis Estimation
#'
#' Estimates time to diagnosis of active TB for each patient using simple pre-defined rules.
#' This time is not directly available in the data but can be got at in most cases indirectly.
#'
#' The rules for estimating time to diagnosis are (in priority order):
#' \enumerate{
#'  \item {If culture +ve then \code{min}(report time, start treatment)
#'        [previously, take report time only.]}
#'  \item {What is the test used ('means') for the final diagnosis? Use corresponding times. Latest date of those for the means of diagnosis used.}
#'  \item If diagnosis date is too late (130 days after first test) then use \code{DateDiagCon}
#'  \item If diagnosis date still \code{NA} then use  \code{DateDiagCon}.
#'  \item {Could use \code{min(DateDiagCon, lasttestDate)} instead?}
#'  }
#'
#' @param data Individual patient records from IDEA dataset
#'
#' @return Original data with the \code{start.to.diag} column appended
#' @export
#'
estimateTimeToDiagnosis <- function(data){


  for (i in 1:nrow(data)){

    if (data$TBcult[i] == "POSITIVE"){

      data$start.to.diag[i] <- min(data$testDrug_diff[i],
                                   data$start.to.TBcultorig[i],
                                   na.rm = TRUE)

      ##TODO: not sure why changed NAs to 0s here?
      # data$start.to.diag[i] <- ifelse(is.na(data$start.to.diag[i]),
      #                                 0,
      #                                 data$start.to.diag[i])

    }else{
      data$start.to.diag[i] <- getDiagnosisTimeViaMeans(data, i)}
  }

  data$start.to.diag[is.infinite(data$start.to.diag)] <- NA
  data$start.to.diag[data$start.to.diag < 0] <- NA
  data$testDiagCon_diff[data$testDiagCon_diff < 0] <- NA

  tooBigLimit <- 130  #days; ~4 months from expert input

  tooBig <-
    data$start.to.diag > tooBigLimit &
    data$testDiagCon_diff < tooBigLimit &
    !is.na(data$testDiagCon_diff)

  start_to_diag_sub <- tooBig | is.na(data$start.to.diag)

  data$start.to.diag[start_to_diag_sub] <- data$testDiagCon_diff[start_to_diag_sub]

  invisible(data)
}

#' estimate_start_to_diag
#'
#' @param data
#'
#' @return
#' @export
#'
estimate_start_to_diag <- function(data) {

  res <- estimateTimeToDiagnosis(data)
  return(res$start.to.diag)
}

#' is.TestUsedAsMeansOfDiag
#'
#' @param data
#' @param i
#' @param testKeyword
#'
#' @return
#' @export
#'
is.TestUsedAsMeansOfDiag <- function(data,
                                     i,
                                     testKeyword){

  # fields to search
  datai <- data[i, c("MeansFinDiag1",
                     "MeansFinDiagother.multiple",
                     "Diagnostic_mean")]

  out <- any(sapply(datai,
             function(x)
               grepl(testKeyword, x, ignore.case = TRUE)))
  return(out)
}


#' getDiagnosisTimeViaMeans
#'
#' @param data
#' @param i
#'
#' @return
#' @export
#'
getDiagnosisTimeViaMeans <- function(data,
                                     i){

  diagtime  <- NA
  diag_list <- list("Smear" = "start.to.Smear",
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

  for (test in names(diag_list)){

    test_used <- is.TestUsedAsMeansOfDiag(data, i, test)

    if (test_used) {

      diagtime <- c(diagtime,
                    data[i, diag_list[[test]]])
    }
  }

  return(max(diagtime, na.rm = TRUE))
}
