#' Sum the cost of all specified treatments for each patient
#'
#' \code{calcPatientCostofTests} takes Access database and updates one of the tables
#'
#' @param testdata dataset test results
#' @return single value per patient
#'


calcPatientCostofTests <- function(testdata){
    ## call: calcPatientCostofTests(data[,c("TBcult","Smear")])

#     testdata <- data.frame(TBcult=c("Not tested","POSITIVE","POSITIVE"),
#                            Smear=c("POSITIVE","Not tested","NEGATIVE"))

    testcosts <- c(TBcult=22.29, #10 #London RC of P of. Tuberculosis: clinical diagnosis and management of tuberculosis, and measures for its prevention and control. 2006;(November 2010):1–308.
                   Smear=1.56, #HTA Systematic review, meta-analysis and economic modelling of molecular diagnostic tests for antibiotic resistance in TB.
                   IGRA=56.24, #Pareek M, Bond M, Shorey J, Seneviratne S, Guy M, White P, et al. Community-based evaluation of immigrant tuberculosis screening using interferon   release assays and tuberculin skin testing: observational study and economic analysis. Thorax. 2012;230–9.
                   TSTres=34.22, ##NCCCC, Tuberculosis: appendices. London: Royal College of Physicians, 2006; National Collaborating Centre for Chronic Conditions. TB Partial Update: Cost-effectiveness analysis of interferon gamma release assay (IGRA) testing for latent tuberculosis. London: NICE, 2010.
                   TSPOT=82, #Local lab cost (IDEA health econ analysis plan)
                   QFN=45, #52 #
                   CXR=28, #NCCCC, Tuberculosis: appendices. London: Royal College of Physicians, 2006;National Collaborating Centre for Chronic Conditions. TB Partial Update: Cost-effectiveness analysis of interferon gamma release assay (IGRA) testing for latent tuberculosis. London: NICE, 2010.
                   CSF=0, #
                   BAL=23.24, #St Marys R & D office
                   HistBiop=25, #St Marys R & D office
                   NeedleAsp=90.21, #St Marys R & D office
                   PCR=202.45, #St Marys R & D office
                   CT=300, #St Marys R & D office
                   MRI=375, #St Marys R & D office
                   PET=713) #St Marys R & D office

    whichTests <- function(x,testdata) colnames(testdata)[x=="POSITIVE" | x=="NEGATIVE" | x=="INDETERMINATE"]
    apply(testdata, 1, function(x) sum(testcosts[whichTests(x,testdata)], na.rm=TRUE))
}

