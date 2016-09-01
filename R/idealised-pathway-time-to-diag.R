
#' Calculate NICE Idealised Time to Diagnosis
#'
#' Uses the pathway frequencies for the NICE idealised pathway (chest X-ray, culture, smear, clinical judgement)
#' and whether the observed Dosanjh final diagnosis is what would be expected from the idealised pathway.
#'
#' For the pathway derived from the data only, this is like a `double-imputation' because we impute missing test results to match-up with the NICE guidelines
#' and then we use imputed times to diagnosis.
#'
#' For when we use predefined values we avoid some of the messiness of the time to test data. At the moment we just use some hard-coded expected values.
#' Future work could use some statistic from the multiple imputations, like the median or mode.
#' We could also (hardcode) data derived averaged instead.
#'
#' @param data IDEA study data set (perhaps bootstrap sampled or imputed).
#' @param use_times_in_data Use the raw times in the IDEA data or the idealised times.
#' @param return Return a list of matrix.
#'
#' @return Times to diagnosis statistics
#' @export
calcIdealisedTimetoDiag <- function(data, BAL=FALSE, Dosanjh3TB=FALSE, use_times_in_data=FALSE, return="list"){

  require(reshape2)
  require(pastecs)

  pathfreq <- data

  if(Dosanjh3TB){
    DosTB <- c(1,2,3); DosnonTB <- 4
  }else{
    DosTB <- c(1,2); DosnonTB <- c(3,4)}

  pathfreq <- within(pathfreq, {
    activeTB <- DosanjhGrouped%in%DosTB
    nonTB <- DosanjhGrouped%in%DosnonTB
  })

  load("../../../output_data/NICEguidelines.RData")
  pathfreq <- merge(pathfreq, NICEguidelines)

  pathfreq$start.to.clinfeatures <- 0

  ## could include averages from total data set instead
  if(!use_times_in_data){
    pathfreq$start.to.CXR <- 1
    pathfreq$start.to.Smear <- 2
    pathfreq$start.to.TBcultorig <- 42
  }

  NICEstart.to.diag <- rep(NA, nrow(pathfreq))

  for (i in 1:nrow(pathfreq)){
    time <- pathfreq[i, as.character(pathfreq$NICEmeansofDiagnosis[i])]
    if(!is.null(time)){
      NICEstart.to.diag[i] <- time
    }
  }

  pathfreq$NICEstart.to.diag <- NICEstart.to.diag

  pathfreq$mismatch <- pathfreq$NICEnonTB!=pathfreq$nonTB

  idealised.timetodiag.total <- sum(pathfreq$NICEstart.to.diag[!pathfreq$mismatch])
  Nidealised <- sum(!pathfreq$mismatch)

  propMatch <- Nidealised/nrow(data)

  mismatches.timetodiag.total <- sum(pathfreq$start.to.diag[pathfreq$mismatch])
  mismatches.timetodiag.total.noNA <- sum(pathfreq$start.to.diag[pathfreq$mismatch], na.rm = T)
  num.noNA <- length(!is.na(pathfreq$start.to.diag[pathfreq$mismatch]))

  if(return=="list"){

    return(list(
      # pathfreq=pathfreq,
      counts=c(total=nrow(data), NICEagree=Nidealised, NICEdisagree=(nrow(data)-Nidealised), propMatch=propMatch),
      times=list(observed=c(totalcost=sum(data$start.to.diag), perperson=mean(data$start.to.diag, na.rm = T)),
                 estimated=list(
                   total=c(mismatches=mismatches.timetodiag.total,
                           noMismatches=idealised.timetodiag.total,
                           total=mismatches.timetodiag.total + idealised.timetodiag.total),
                   perperson=c(mismatches=mismatches.timetodiag.total.noNA/num.noNA,
                               noMismatches=idealised.timetodiag.total/Nidealised,
                               total=(mismatches.timetodiag.total.noNA + idealised.timetodiag.total)/(num.noNA+Nidealised)))
      )))
  }else{
    ## for paper table
    return(round(c(
      nrow(data),
      Nidealised,
      (nrow(data)-Nidealised),                                  #mismatched
      propMatch,
      mean(data$start.to.diag, na.rm = T),                      #mean all observed (per patient)
      idealised.timetodiag.total/Nidealised,                    #mean idealised
      mismatches.timetodiag.total.noNA/num.noNA,                #mean mismatched observed
      (mismatches.timetodiag.total.noNA + idealised.timetodiag.total)/(num.noNA+Nidealised)),2)) #mean mismatched observed and idealised
  }
}

