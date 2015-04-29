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
                  "step1Diag","Sx_resolved","treatResponse","Alt_diag")

  pdata <- prepDataForDecTree(pdata, data.names)

  out <- createDecisionTree(pdata, test.DT, diag_obs.DT, diag_est.DT)

  tabs <- writeAllFreqTables(out, dir=dir)

  tabs
}

