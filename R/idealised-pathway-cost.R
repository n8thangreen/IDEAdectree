
#' Calculate NICE idealised test costs
#'
#' Uses the pathway frequencies for the NICE idealised pathway (chest X-ray, culture, smear, clinical judgement)
#' and whether the observed Dosanjh final diagnosis is what would be expected from the idealised pathway.
#' Where it is not (usually when the 4 measurement are insufficient to discriminate) then we use
#' the raw data and just sum the cost of the observed tests. The subjectivity come in when deciding
#' exactly what the expected diagnosis is. In most cases it is clear. We remove patients who do not have all 3 tests
#' so there may be some bias introduced.
#'
#' This is the function used in \code{\link{boot.calcIdealisedTestCosts}}.
#'
#' @param data IDEA data set (perhaps bootstrap sampled)
#'
#' @return Total test cost by Dosanjh category
#'

calcIdealisedTestCosts <- function(data){

  require(reshape2)
  require(pastecs)

  ## distiguish between when clinical features used to rule-in or rule-out TB
  clinfeatures <- rep("none", nrow(data))
  clinfeatures[grepl(x = data$Diagnostic_mean, pattern="Clinical") & data$TBconfirmed=="TRUE"]  <- "TB"
  clinfeatures[grepl(x = data$Diagnostic_mean, pattern="Clinical") & data$TBconfirmed=="FALSE"] <- "other"

  dat <- cbind(ones=1, clinfeatures, data[,c("PatientStudyID","TBcult","Smear","CXR","DosanjhGrouped")])
  pathfreq <- dcast(dat, TBcult+Smear+CXR+clinfeatures~DosanjhGrouped, value.var = "ones", fun.aggregate = length)

  ## had at least Culture, CXR and Smear
  pathfreq <- pathfreq[!apply(pathfreq,1, function(y) "Not taken"%in%y),]

  ## pathway frequencies without the `unusual'/wrongly diagnosed patients
  x.without <- pathfreq
  Nx <- nrow(x.without)

  ## pathways with multiple outcomes
  which4  <- x.without$TBcult=="INDETERMINATE" & x.without$Smear=="POSITIVE" & x.without$CXR=="POSITIVE" & x.without$clinfeatures=="none"
  which21 <- x.without$TBcult=="NEGATIVE" & x.without$Smear=="NEGATIVE" & x.without$CXR=="INDETERMINATE" & x.without$clinfeatures=="none"
  which22 <- x.without$TBcult=="NEGATIVE" & x.without$Smear=="NEGATIVE" & x.without$CXR=="NEGATIVE" & x.without$clinfeatures=="none"
  which23 <- x.without$TBcult=="NEGATIVE" & x.without$Smear=="NEGATIVE" & x.without$CXR=="POSITIVE" & x.without$clinfeatures=="none"

  Nunusual <- c(0,0,0,0)
  Nunusual[4] <- ifelse(sum(which4)>0, pathfreq[which4,"4"], 0)
  Nunusual[2] <- pathfreq[which21,"2"] + pathfreq[which22,"2"] + pathfreq[which23,"2"]

  Nunusual.total <- sum(Nunusual)

  x.without[which4,"4"]  <- 0
  x.without[which21,"2"] <- 0
  x.without[which22,"2"] <- 0
  x.without[which23,"2"] <- 0

  categories <- c("1","2","3","4")
  x.without <- cbind(x.without, sum=rowSums(pathfreq[,categories]))
  x.without <- rbind(x.without, c("","","","",colSums(x.without[,categories])))
  Nx <- Nx + 1

  idealisedTestsCosts <- c(TBcult=22.29, Smear=1.56, CXR=28)
  idealisedTestCost.total  <- sum(idealisedTestsCosts)
  Nidealised.Dosanjh <- as.numeric(x.without[Nx,categories])
  Nidealised <- sum(Nidealised.Dosanjh)

  ## cost for patients with idealised tests and `correct' diagnosis
  idealisedpatientsum.Dosanjh <- Nidealised.Dosanjh * idealisedTestCost.total
  idealisedpatientsum <- sum(idealisedpatientsum.Dosanjh)


  ## + test cost for patients with `incorrect' diagnosis
  unusualcost <- c(0,0,0,0)

  patient4 <- data$TBcult==pathfreq[which4,"TBcult"] &
    data$Smear==pathfreq[which4,"Smear"] &
    data$CXR==pathfreq[which4,"CXR"] &
    data$DosanjhGrouped=="4"

  patient21 <- data$TBcult==pathfreq[which21,"TBcult"] &
    data$Smear==pathfreq[which21,"Smear"] &
    data$CXR==pathfreq[which21,"CXR"] &
    data$Diagnostic_mean=="Imaging" &
    data$DosanjhGrouped=="2"

  patient22 <- data$TBcult==pathfreq[which22,"TBcult"] &
    data$Smear==pathfreq[which22,"Smear"] &
    data$CXR==pathfreq[which22,"CXR"] &
    data$Diagnostic_mean!="Clinical Features" &
    data$DosanjhGrouped=="2"

  patient23 <- data$TBcult==pathfreq[which23,"TBcult"] &
    data$Smear==pathfreq[which23,"Smear"] &
    data$CXR==pathfreq[which23,"CXR"] &
    data$Diagnostic_mean!="Clinical Features" &
    data$DosanjhGrouped=="2"

  listunusualcost <- vector("list",length = 4)

  if(Nunusual[4]>0){
    listunusualcost[[4]] <- data$totalcost[patient4]
    unusualcost[4] <- sum(listunusualcost[[4]])
  }
  if(Nunusual[2]>0){
    listunusualcost[[2]] <- data$totalcost[patient21]
    listunusualcost[[2]] <- c(listunusualcost[[2]], data$totalcost[patient22])
    listunusualcost[[2]] <- c(listunusualcost[[2]], data$totalcost[patient23])

    unusualcost[2] <- sum(listunusualcost[[2]])
  }

  ## create idealised cost sample
  for (i in 1:4){
    listunusualcost[[i]] <- c(rep(idealisedTestCost.total, Nidealised.Dosanjh[[i]]), listunusualcost[[i]])
  }

  ##total MEAN test cost by Dosanjh category
  res <- NULL
  for (i in as.numeric(categories)){

    res[i] <- round((unusualcost[i] + idealisedpatientsum.Dosanjh[i])/(Nunusual[i] + Nidealised.Dosanjh[i]),2)
  }

  res <- list(means=res, stats=sapply(listunusualcost, summary))
  return(res)
}


#' Bootstrap calculate idealised test costs
#'
#'
#'
#' @param data IDEA data set
#' @param stat statistic (Mean, Median, 1st Qu., 3rd Qu.)
#' @param n number of bootstrap samples
#'
#' @return Array of n bootstrap outputs
#'
#' See also \code{\link{calcIdealisedTestCosts}}.
#'

boot.calcIdealisedTestCosts <- function(data, stat="Mean", n=100){

  out <- NULL
  for (i in 1:n){
    data.boot <- data[sample(1:nrow(data), replace = TRUE, size = nrow(data)),]

    out <- rbind(out, calcIdealisedTestCosts(data.boot)$stats[stat,])
    # means <- rbind(means, calcIdealisedTestCosts(data.boot)$means)
  }
  return(out)
}



