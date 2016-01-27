
#' Sum the Cost of All Specified Treatments for Each Patient
#'
#' \code{calcPatientCostofTests} sums the total cost of diagnostic tests in the pathway.
#'
#' If a BAL is taken we include the cost of the proceduce and the additional cost of a culture and smear if used for these.
#'
#' @param data IDEA study data set
#' @param COSTS Named vector each test/procedure and unit costs.
#' @param x3 Should the smear and culture be multiplied by 3? It is standard practice to perform 3 of each with a sputum sample.
#'
#' @return single value per patient

calcPatientCostofTests <- function(data, COSTS, x3 = TRUE){

  freq <- ifelse(x3, 3, 1)
  if(is(COSTS, "list")){COSTS <- COSTS[["mean"]]}

  namesNumEachTest <- grep(names(data), pattern = "Num_",value = TRUE)

  testcosts <- c(TBcult = freq*COSTS$TBcult,
                 Smear = freq*COSTS$Smear,
                 IGRA=COSTS$IGRA,
                 TSTres=COSTS$TSTres,
                 TSPOT=COSTS$TSPOT,
                 QFN=COSTS$QFN,
                 CXR=COSTS$CXR,
                 CSF=COSTS$CSF,
                 BAL=COSTS$BAL,
                 HistBiop=COSTS$HistBiop,
                 NeedleAsp=COSTS$NeedleAsp,
                 PCR=COSTS$PCR,
                 CT=COSTS$CT,
                 MRI=COSTS$MRI,
                 PET=COSTS$PET,
                 TBcultSmearBAL = COSTS$TBcult + COSTS$Smear)#,
                 # Risk_factors=0)  #why did I include this??

  ## if there's no test count data
  if(length(namesNumEachTest)==0){
    data <- as.matrix(data) #need to do it this way because factor
    data[data%in%c("POSITIVE", "NEGATIVE")] <- 1
    data[data!=1 | is.na(data)] <- 0
    class(data) <- "numeric"
    testcosts[colnames(data)]%*%t(data)
  }else{
    ## cross-product to account for multiple same tests per patient
    c(testcosts%*%t(data[,paste("Num_",names(testcosts),sep="")]))
  }
}


#' Sample from Standard Distributions
#'
#' @param param.distns
#'
#' @return vector of sample points
#' @examples
#'
#' param.distns <- list(TBcult=list(distn="gamma",
#' params=c(mean=22.29, sd=2.23)),
#' Smear=list(distn="gamma",
#'            params=c(mean=7, sd=0.68)),
#' IGRA=list(distn="unif",
#'           params=c(min=24, max=100)),
#' TSTres=list(distn="unif",
#'             params=c(min=8, max=32)),
#' TSPOT=list(distn="unif",
#'            params=c(min=45, max=99)),
#' QFN=list(distn="unif",
#'          params=c(min=36.8, max=84)),
#' CXR=list(distn="unif",
#'          params=c(min=26, max=41)),
#' CSF=list(distn="none",
#'          params=c(mean=0)),
#' BAL=list(distn="none",
#'          params=c(mean=612+23.24)),
#' HistBiop=list(distn="none",
#'               params=c(mean=25)),
#' NeedleAsp=list(distn="none",
#'                params=c(mean=90.21)),
#' PCR=list(distn="none",
#'          params=c(mean=202.45)),
#' CT=list(distn="none",
#'         params=c(mean=300)),
#' MRI=list(distn="none",
#'          params=c(mean=375)),
#' PET=list(distn="none",
#'          params=c(mean=713)))

sample.distributions <- function(param.distns){

  out <- data.frame(matrix(NA, nrow = 1, ncol = length(COST.distns)))
  for (i in 1:length(param.distns)){

    out[i] <- switch(param.distns[[i]]$distn,
                     gamma = rgamma(1, shape = MoM.gamma(mean=param.distns[[i]]$params["mean"],
                                                         var=param.distns[[i]]$params["sd"]^2)$shape,
                                       scale = MoM.gamma(mean=param.distns[[i]]$params["mean"],
                                                         var=param.distns[[i]]$params["sd"]^2)$scale),
                     unif = runif(1, param.distns[[i]]$params["min"],
                                     param.distns[[i]]$params["max"]),
                     none = param.distns[[i]]$params["mean"])
  }

  names(out) <- names(param.distns)

  return(out)
}

