

#' calc.currentDiagTimeAndCost.bootsample
#'
#' bootstrap patient data for mean and standard deviation estimates
#' for current in-practice routine pathway (without rule-out test)
#'
#' @param data
#'
#' @return bootstrapped dataset
#' @export
#'

calc.currentDiagTimeAndCost.bootsample <- function(data){

  bootIDs <- sample(1:nrow(data), replace = TRUE)
  data[bootIDs,]
}


#' calc.enhancedDiagTimeAndCost.bootsample
#'
#' bootstrap patient data for mean and standard deviation estimates
#' for enhanced pathway (with rule-out test) by representing as mixture model
#'
#' @param data
#'
#' @return bootstrapped dataset
#' @export
#'

calc.enhancedDiagTimeAndCost.bootsample <- function(data){

  ECDF.TB <- ecdf(data$riskfacScore[data$DosanjhGrouped%in%c(1,2)])
  ECDF.nonTB <- ecdf(data$riskfacScore[data$DosanjhGrouped==4])

  prop_highrisk <- 0.4
  gamma.TB <- 1-ECDF.TB(prop_highrisk)
  gamma.nonTB <- 1-ECDF.nonTB(prop_highrisk)

  spec  <- 0.9
  sens  <- 0.9
  FNtime <- 42

  testcost <- 200

  # randomly select patients for each pathway
  TBids <- which(data$DosanjhGrouped%in%c(1,2,3))
  nonTBids <- which(data$DosanjhGrouped%in%c(4))

  ## for TB
  whosleft <- TBids
  nleft <- length(whosleft)

  clinjudge <- sample(x=whosleft, size=gamma.TB*nleft, replace=FALSE)
  whosleft  <- whosleft[!whosleft%in%clinjudge]
  nleft <- length(whosleft)

  testpos <- sample(x=whosleft, size=nleft*sens, replace=FALSE)
  data$start.to.diag[testpos] <- data$start.to.diag[testpos] + 1
  data$totalcost[testpos] <- data$totalcost[testpos] + testcost
  whosleft  <- whosleft[!whosleft%in%testpos]
  nleft <- length(whosleft)

  delay <- whosleft
  data$start.to.diag[delay] <- data$start.to.diag[delay] + 1 + FNtime
  data$totalcost[delay] <- data$totalcost[delay] + testcost


  ## for nonTB
  whosleft <- nonTBids
  nleft <- length(whosleft)

  clinjudge <- sample(x=whosleft, size=gamma.nonTB*nleft, replace=FALSE)
  whosleft  <- whosleft[!whosleft%in%clinjudge]
  nleft <- length(whosleft)

  testpos <- sample(x=whosleft, size=nleft*(1-spec), replace=FALSE)
  data$start.to.diag[testpos] <- data$start.to.diag[testpos] + 1
  data$totalcost[testpos] <- data$totalcost[testpos] + testcost
  whosleft  <- whosleft[!whosleft%in%testpos]
  nleft <- length(whosleft)

  delay <- whosleft
  data$start.to.diag[delay] <- 1
  data$totalcost[delay] <- testcost

  data
}


#' make.tableDiagCost_bootmeanse
#'
#' bootstrap patient data for mean and standard deviation estimates
#'
#' @param data
#' @param sampleFUN
#'
#' @return matrix of means and sd
#' @export
#'

make.tableDiagCost_bootmeanse <- function(data, sampleFUN=calc.enhancedDiagTimeAndCost.bootsample){

  iterations <- 1000
  table <- array(NA, dim=c(4,13,iterations))
  for (i in 1:iterations){
    table[,,i] <- maketable.Dosanjh_TimeCost(sampleFUN(data))
  }

  meanvals <- apply(table, c(1,2), function(x) round(mean(as.numeric(x)),2))
  sdvals <- apply(table, c(1,2), function(x) round(sd(as.numeric(x)),2))
  matrix(paste(meanvals, " (",sdvals,")", sep=""), nrow = 4)
}






