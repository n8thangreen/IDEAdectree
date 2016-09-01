
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
#' @return Single value per patient
#' @export
calcPatientCostofTests <- function(data, COSTS, x3 = TRUE){

  require(assertive)
  x3 <- coerce_to(use_first(x3), "logical")

  freq <- ifelse(x3, 3, 1)
  if(is(COSTS, "list")){COSTS <- COSTS[["mean"]]}

  namesNumEachTest <- grep(names(data), pattern = "Num_",value = TRUE)

  testcosts <- c(TBcult = freq*COSTS$TBcult,
                 Smear = freq*COSTS$Smear,
                 IGRA = (COSTS$TSPOT+COSTS$QFN)/2,  #COSTS$IGRA,
                 TSTres = COSTS$TSTres,
                 TSPOT =COSTS$TSPOT,
                 QFN = COSTS$QFN,
                 CXR = COSTS$CXR,
                 CSF = COSTS$CSF,
                 BAL = COSTS$BAL,
                 HistBiop = COSTS$HistBiop,
                 NeedleAsp = COSTS$NeedleAsp,
                 PCR = COSTS$PCR,
                 CT = COSTS$CT,
                 MRI = COSTS$MRI,
                 PET = COSTS$PET,
                 EBUS = COSTS$EBUS,
                 TBcultSmearBAL = COSTS$TBcult + COSTS$Smear)#,
                 # Risk_factors=0)  #why did I include this??

  ## if there's no test count data
  if(length(namesNumEachTest)==0){

    data <- as.matrix(data) #need to do it this way because factor
    data[data%in%c("POSITIVE", "NEGATIVE")] <- 1L  #i.e. taken
    data[data!=1 | is.na(data)] <- 0L
    class(data) <- "numeric"
    testcosts[colnames(data)]%*%t(data)
  }else{
    ## cross-product to account for multiple same tests per patient
    Num_names <- paste("Num_", names(testcosts), sep="")
    return(c(testcosts%*%t(data[ ,Num_names])))
  }
}


#' Sample from Standard Distributions
#'
#' Supply a list of defined distributions and one sample realisation is taken of each.
#'
#' @param param.distns List of distribution names and their respective parameter values
#'
#' @return vector of sample points
#' @examples
#' @export
#' ## input used in IDEA HTA C-E analysis
#'
#'  param.distns <- list(TBcult=list(distn="gamma",
#'  params=c(mean=22.29, sd=2.23)),
#'  Smear=list(distn="gamma",
#'             params=c(mean=7, sd=0.68)),
#'  IGRA=list(distn="unif",
#'            params=c(min=15, max=58)),
#'  TSTres=list(distn="unif",
#'              params=c(min=8, max=36)),
#'  TSPOT=list(distn="unif",
#'             params=c(min=50, max=106)),
#'  QFN=list(distn="unif",
#'           params=c(min=26, max=78)),
#'  CXR=list(distn="unif",
#'           params=c(min=23, max=43)),
#'  CSF=list(distn="none",
#'           params=c(mean=0)),
#'  BAL=list(distn="unif",
#'           params=c(min=306, max=1224)),
#'  HistBiop=list(distn="unif",
#'                params=c(min=12.5, max=50)),
#'  NeedleAsp=list(distn="unif",
#'                 params=c(min=45.1, max=180.42)),
#'  PCR=list(distn="unif",
#'           params=c(min=101.2, max=404.9)),
#'  CT=list(distn="none",
#'          params=c(mean=300)),
#'  MRI=list(distn="none",
#'           params=c(mean=375)),
#'  PET=list(distn="none",
#'           params=c(mean=713)),
#'  EBUS=list(distn="none",
#'           params=c(mean=2634)))

sample.distributions <- function(param.distns){

  require(triangle)

  out <- data.frame(matrix(NA, nrow = 1, ncol = length(param.distns)))
  for (i in 1:length(param.distns)){

    out[i] <- switch(param.distns[[i]]$distn,
                     gamma = {
                              mom <- MoM.gamma(mean=param.distns[[i]]$params["mean"], var=param.distns[[i]]$params["sd"]^2)
                              gamma = rgamma(1, shape = mom$shape, scale = mom$scale)
                     },
                     unif = runif(1, param.distns[[i]]$params["min"],
                                  param.distns[[i]]$params["max"]),
                     triangle = rtriangle(1, param.distns[[i]]$params["min"],
                                          param.distns[[i]]$params["max"]),
                     none = param.distns[[i]]$params["mean"])
  }

  names(out) <- names(param.distns)

  return(out)
}

