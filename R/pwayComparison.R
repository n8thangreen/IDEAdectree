#' Compare data and idealised care pathways.
#'
#' \code{pwayComparison} returns joined tables of inputs and predicted outcomes
#'
#' @param data Individual patient records
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

pwayComparison <- function(pdata, dir="C:/Users/Nathan/Dropbox/TB/IDEA/output_data/from IDEAdectree"){

  #if(!interactive)?

  data.names <- c("PatientStudyID","TBconfirmed","DosanjhGrouped",
                  "TBcult","Smear","TSTres","CXR","IGRA",
                  "CXR","QFN","TSPOT","CT","PET","MRI","PCR",
                  "step1Diag","Sx_resolved","treatResponse","Alt_diag","timeToDiag")

  pdata <- prepDataForDecTree(pdata, data.names)

  out <- createDecisionTree(pdata, test.DT, diag_obs.DT, diag_est.DT)

  tabs <- writeAllFreqTables(out, dir=dir)

  tabs
}


#' A simplified version of comparison of care pathways.
#'
#' \code{pwayComparison.simple} returns one pathway aggregated table and median time to diagnoses
#'
#' @param data Individual patient records
#' @param test.DT (idealised) decision tree object returned by \code{createDTobj()}
#' @param form formula identifier for grouping by final Dosanjh category (default), observed presumptive treatment decision (step1Diag) or pooled (all)
#' @return list length 2

pwayComparison.simple <- function(data, test.DT, form="Dosanjh"){

  data.names <- c("PatientStudyID","TBconfirmed","DosanjhGrouped",
                  "TBcult","Smear","TSTres","CXR","IGRA",
                  "CXR","QFN","TSPOT","CT","PET","MRI","PCR",
                  "step1Diag","Sx_resolved","treatResponse","Alt_diag","timeToDiag")

  data <- prepDataForDecTree(data, data.names)

  formula <- switch(form,
                    all = as.formula(paste(paste(names(test.DT$input), collapse="+"),"~ .")),
                    Dosanjh = as.formula(paste(paste(names(test.DT$input), collapse="+"),"~ DosanjhGrouped")),
                    step1Diag = as.formula(paste(paste(names(test.DT$input), collapse="+"),"~ step1Diag")))

  data$Freq <- 1
  dt <- reshape::cast(data, formula, fun=sum, value="Freq")
  eachCost <- calcPatientCostofTests(dt[,names(test.DT$input)])
  eachCost[length(eachCost)] <- sum(eachCost)
  dt <- data.frame(dt, Costs = eachCost)
  dt
  # write.csv(dt, file=paste(dir,"/testsonly.csv", sep=""))

  dt.medianTimeToDiag <- reshape::cast(data, formula, fun=median, value="timeToDiag", na.rm=T)

  # write.csv(dt.medianTimeToDiag, file=paste(dir,"/median-time-to-diagnosis.csv", sep=""))
  ##TODO##
  #doesn't account for missing (censored?) times tho

  list(dt=dt, dt.medianTimeToDiag=dt.medianTimeToDiag)
}


