#' Rule-based time to diagnosis estimation
#'
#' \code{estimateTimeToDiagnosis} returns estimates of time to diagnosis for each patient
#'
#' @param data Individual patient records
#' @return timeToDiag

estimateTimeToDiagnosis <- function(data){
  ## to estimate the length of time to final (confirmed) diagnosis
  ##RULES

  ## time of diagnosis field is upper limit
  TimeToDiag <- data$testDiagCon_diff

  for (i in 1:nrow(data)){

    if (data$TBconfirmed[i]==TRUE & data$TBDrug_diff[i]>0 & !is.na(data$TBDrug_diff[i])){
      ## assume that all TBconfirmed to eventually go on treatment!

      ## to end of treatment
      t2d <- data$testDrug_diff[i] + data$TBDrug_diff[i]

      ## first date decision to be kept on treatment e.g. 2 months after starting (assuming the response to treatment)
      # difftime(add(TBDrugStart.min,2), testDate.min)

      ## first date started on treatment?
      # difftime(TBDrugStart.min, testDate.min)
      # date of followup if on treatment afterwards ie enddate treatment>followup date
      # date of (all?) symptoms resolving given on treatment
    }

    if (data$TBconfirmed[i]==FALSE){
      ## put on treatment only
      ## last date that time taken off treatment | started treatment in the first place
      if(data$TBDrug_diff[i]>0 & !is.na(data$TBDrug_diff[i])){
        t2d <- data$testDrug_diff[i] + data$TBDrug_diff[i]
      }else{
        ## not dependent on treatment
        #       time of alternative diagnosis
        #       time of start treatment for alternative diagnosis

        t2d <- NA
      }

    }

    if (!is.na(TimeToDiag[i]) & !is.na(t2d)){
      TimeToDiag[i] <- min(as.numeric(TimeToDiag[i]), t2d, na.rm = T)}
  }

  as.numeric(TimeToDiag)
}



