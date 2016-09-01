
#' Calculate NICE Idealised Test Costs
#'
#' Uses the pathway frequencies for the NICE idealised pathway (chest X-ray, culture, smear, clinical judgement)
#' and whether the observed Dosanjh final diagnosis is what would be expected from the idealised pathway.
#' Assume that three sputum samples are taken so at least 3 culture and smear tests.
#' Inclusion of BAL would further include this procedure cost and the additional culture and smear tests.
#' Where it is not (usually when the 4 measurement are insufficient to discriminate) then we use
#' the raw data and just sum the cost of the observed tests. The subjectivity come in when deciding
#' exactly what the expected diagnosis is. In most cases it is clear.
#' We could remove patients who do not have all 3 tests so there may be some bias introduced,
#' but can use an imputed completed sample (prefered).
#'
#' This is the function used in \code{\link{boot.calcIdealisedTestCosts}}.
#'
#' @param data IDEA data set (perhaps bootstrap sampled or imputed)
#' @param return Return a list of matrix.
#'
#' @return Total test cost statistics

calcIdealisedTestCosts <- function(data, BAL=FALSE, Dosanjh3TB=FALSE, return="matrix"){

  require(reshape2)
  require(pastecs)

  pathfreq <- data[,c("TBcult","Smear","CXR","clinfeatures","DosanjhGrouped","totalcost")]

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

  pathfreq$mismatch <- pathfreq$NICEnonTB!=pathfreq$nonTB

  ## three sputum samples taken. BAL procedure+test
  if(BAL){
    idealisedTestsCosts <- c(TBcult=(612+23.24) + (22.29*4), Smear=1.56*4, CXR=28)
  }else{
    idealisedTestsCosts <- c(TBcult=22.29*3, Smear=1.56*3, CXR=28)}

  idealisedTestCost.total <- sum(idealisedTestsCosts)
  Nidealised <- sum(!pathfreq$mismatch)

  totalTestCost.noMismatches <- Nidealised*idealisedTestCost.total  #for original complete data
  ##TODO## for incomplete data

  propMatch <- Nidealised/nrow(data)

  totalTestCost.mismatches <- sum(pathfreq$totalcost[pathfreq$mismatch])

  if(return=="list"){
  return(list(
    # pathfreq=pathfreq,
    counts=c(total=nrow(data), NICEagree=Nidealised, NICEdisagree=(nrow(data)-Nidealised), propMatch=propMatch),
    costs=list(observed=c(totalcost=sum(data$totalcost), perperson=sum(data$totalcost)/nrow(data)),
               estimated=list(
                 total=c(mismatches=totalTestCost.mismatches,
                         noMismatches=totalTestCost.noMismatches,
                         total=totalTestCost.mismatches + totalTestCost.noMismatches),
                 perperson=c(mismatches=totalTestCost.mismatches/(nrow(data)-Nidealised),
                             noMismatches=totalTestCost.noMismatches/Nidealised,
                             total=(totalTestCost.mismatches + totalTestCost.noMismatches)/nrow(data)))
    )))
  }else{
  ## for paper table
  return(round(c(
    sum(data$totalcost),                                                   #total all observed
    sum(data$totalcost)/nrow(data),                                        #mean all observed (per patient)
    totalTestCost.noMismatches,                                            #total idealised
    totalTestCost.noMismatches/Nidealised,                                 #mean idealised
    totalTestCost.mismatches,                                              #total mismatched observed
    totalTestCost.mismatches/(nrow(data)-Nidealised),                      #mean mismatched observed
    totalTestCost.mismatches + totalTestCost.noMismatches,                 #total mismatched observed and idealised
    (totalTestCost.mismatches + totalTestCost.noMismatches)/nrow(data)),0))#mean mismatched observed and idealised
  }
}


#' Bootstrap calculate idealised test costs [DEPRECATED]
#'
#' @param data IDEA data set
#' @param stat Statistic (Mean, Median, 1st Qu., 3rd Qu.)
#' @param n Number of bootstrap samples
#'
#' @return Array of \code{n} bootstrap outputs
#'
#' @seealso \code{\link{calcIdealisedTestCosts}}.
#' @export

boot.calcIdealisedTestCosts <- function(data, stat="Mean", n=100){

  out <- NULL
  for (i in 1:n){
    data.boot <- bootstrap.data(data)
    out <- rbind(out, calcIdealisedTestCosts(data.boot)$stats[stat,])
    # means <- rbind(means, calcIdealisedTestCosts(data.boot)$means)
  }
  return(out)
}



