
#' Calculate NICE idealised test costs
#'
#' Uses the pathway frequencies for the NICE idealised pathway (chest X-ray, culture, smear, clinical judgement)
#' and whether the observed Dosanjh final diagnosis is what would be expected from the idealised pathway.
#' Assume that three sputum samples are taken so at at least 3 culture and smear tests.
#' Inclusion of BAL would durther include this procedure cost and the additional culture and smear tests.
#' Where it is not (usually when the 4 measurement are insufficient to discriminate) then we use
#' the raw data and just sum the cost of the observed tests. The subjectivity come in when deciding
#' exactly what the expected diagnosis is. In most cases it is clear.
#' We remove patients who do not have all 3 tests so there may be some bias introduced, but can use an imputed completed sample.
#'
#' This is the function used in \code{\link{boot.calcIdealisedTestCosts}}.
#'
#' @param data IDEA data set (perhaps bootstrap sampled or imputed)
#'
#' @return Total test cost by Dosanjh category

calcIdealisedTestCosts <- function(data){

  require(reshape2)
  require(pastecs)

#   dat <- cbind(ones=1, data[,c("PatientStudyID","TBcult","Smear","CXR","clinfeatures","DosanjhGrouped")])
#
#   pathfreq <- dcast(dat, formula = TBcult+Smear+CXR+clinfeatures~DosanjhGrouped, value.var = "ones", fun.aggregate = length)
#
#   pathfreq <- within(pathfreq, {
#     activeTB <- `1`+`2`
#     nonTB <- `4`
#     rm(`1`,`2`,`3`,`4`)
#   })
#
#   ## remove incomplete data i.e. at least Culture, CXR and Smear
#   pathfreq <- pathfreq[!apply(pathfreq,1, function(y) "Not taken"%in%y),]

  pathfreq <- data[,c("PatientStudyID","TBcult","Smear","CXR","clinfeatures","DosanjhGrouped")]

  pathfreq <- within(pathfreq, {
    activeTB <- DosanjhGrouped%in%c(1,2)
    nonTB <- DosanjhGrouped%in%c(3,4)
  })

  ## pathway frequencies without wrongly diagnosed patients
  pathfreq.noMismatches <- pathfreq
  Npaths <- nrow(pathfreq.noMismatches)

  load("C:/Users/ngreen1/Dropbox/TB/IDEA/output_data/NICEguidelines.RData")
  pathfreq <- merge(pathfreq, NICEguidelines)

  ## mismatches
#   mismatch.nonTB <- pathfreq$NICEnonTB==FALSE & pathfreq$nonTB!=0
#   mismatch.activeTB <- pathfreq$NICEactiveTB==FALSE & pathfreq$activeTB!=0

  mismatch <- pathfreq$NICEnonTB!=pathfreq$nonTB

  Nmismatch <- c(activeTB=0, nonTB=0)
  Nmismatch["activeTB"] <- sum(pathfreq[mismatch.activeTB, "activeTB"])
  Nmismatch["nonTB"] <- sum(pathfreq[mismatch.nonTB, "nonTB"])
  Nmismatch.total <- sum(Nmismatch)

  pathfreq.noMismatches[mismatch.activeTB,"activeTB"] <- 0
  pathfreq.noMismatches[mismatch.nonTB,"nonTB"] <- 0

  NactiveTB.noMismatches <- sum(pathfreq.noMismatches$activeTB)
  NnonTB.noMismatches <- sum(pathfreq.noMismatches$nonTB)

  idealisedTestsCosts <- c(TBcult=22.29*3, Smear=1.56*3, CXR=28)  #three sputum samples taken
  # idealisedTestsCosts <- c(TBcult=23.24+22.29*4, Smear=1.56*4, CXR=28)  #plus BAL

  idealisedTestCost.total <- sum(idealisedTestsCosts)
  Nidealised <- NactiveTB.noMismatches + NnonTB.noMismatches

  totalTestCost.noMismatches <- Nidealised*idealisedTestCost.total  #for complete data
  ##TODO## for incomplete data

  Nidealised/nrow(data)






#   ## + test cost for patients with `incorrect' diagnosis
#   cost.mismatched <- c(activeTB=0, nonTB=0)
#   mismatch.nonTB
#   listunusualcost <- vector("list",length = 4)
#
#   if(Nunusual[4]>0){
#     listunusualcost[[4]] <- data$totalcost[patient4]
#     unusualcost[4] <- sum(listunusualcost[[4]])
#   }
#   if(Nunusual[2]>0){
#     listunusualcost[[2]] <- data$totalcost[patient21]
#     listunusualcost[[2]] <- c(listunusualcost[[2]], data$totalcost[patient22])
#     listunusualcost[[2]] <- c(listunusualcost[[2]], data$totalcost[patient23])
#
#     unusualcost[2] <- sum(listunusualcost[[2]])
#   }
#
#   ## create idealised cost sample
#   for (i in 1:4){
#     listunusualcost[[i]] <- c(rep(idealisedTestCost.total, Nidealised.Dosanjh[[i]]), listunusualcost[[i]])
#   }
#
#   ## total MEAN test cost by Dosanjh category
#   res <- NULL
#   for (i in as.numeric(categories)){
#
#     res[i] <- round((unusualcost[i] + idealisedpatientsum.Dosanjh[i])/(Nunusual[i] + Nidealised.Dosanjh[i]),2)
#   }
#
#   res <- list(means=res, stats=sapply(listunusualcost, summary))
  return(res)
}


#' Bootstrap calculate idealised test costs
#'
#' @param data IDEA data set
#' @param stat Statistic (Mean, Median, 1st Qu., 3rd Qu.)
#' @param n Number of bootstrap samples
#'
#' @return Array of \code{n} bootstrap outputs
#'
#' @seealso \code{\link{calcIdealisedTestCosts}}.
#'

boot.calcIdealisedTestCosts <- function(data, stat="Mean", n=100){

  out <- NULL
  for (i in 1:n){
    data.boot <- bootstrap.data(data)
    out <- rbind(out, calcIdealisedTestCosts(data.boot)$stats[stat,])
    # means <- rbind(means, calcIdealisedTestCosts(data.boot)$means)
  }
  return(out)
}



