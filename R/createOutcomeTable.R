
#' create.indivpathmatrix
#'
#' \code{create.indivpathmatrix} create a table of each distinct pathway through the decision tree
#'
#' @param data
#' @return pathdata
#'

create.indivpathmatrix <- function(data){

  test.names <- c("IGRA", "Smear", "TBcult",# "QFN", "TSPOT", "CSF", "BAL",
                  "HistBiop", "PCR", "imaging")# "PET", "CXR", "NeedleAsp", "CT", "MRI")
  pathdata <- data[,test.names]
  names(pathdata) <- paste(names(pathdata), "bin", sep="_")
  pathdata <- as.matrix(pathdata)
  pathdata[pathdata!="Not taken"] <- 1
  pathdata[pathdata=="Not taken"] <- 0
  pathdata
}


#' get.uniqueDataTestCombinations
#'
#' \code{get.uniqueDataTestCombinations}
#'
#' @param pathdata
#' @return pathdata
#'

get.uniqueDataTestCombinations <- function(pathdata){

  pathdata <- as.data.frame(pathdata[!duplicated(pathdata),])
  pathdata <- pathdata[do.call(order, pathdata),]
  pathdata <- data.frame(pathID=1:nrow(pathdata), pathdata)
  pathdata
}


#' calc.statOutcomeByCovariate
#'
#' \code{calc.statOutcomeByCovariate}
#'
#' @param field
#' @param value
#' @param stat
#' @param outcome
#' @return data
#'

calc.statOutcomeByCovariate <- function(field, value, stat, outcome){
  form <- as.formula(paste(outcome, "~ pathID"))
  aggregate(form, data=data.joined[data.joined[,field]%in%value,], stat, na.rm=T)
}


#' get.meansdOutcomeByCovariate
#'
#' \code{get.meansdOutcomeByCovariate}
#'
#' @param field
#' @param value
#' @param cutoff
#' @param outcome
#' @return data

get.meansdOutcomeByCovariate <- function(field, value, cutoff=8, outcome){

  xbar <- calc.statOutcomeByCovariate(field, value, mean, outcome)
  stdev <- calc.statOutcomeByCovariate(field, value, sd, outcome)
  x <- merge(xbar, stdev, by="pathID", all=TRUE)
  names(x)[2] <- paste(value, "mean", sep="_")
  names(x)[3] <- paste(value, "sd", sep="_")

  quants <- as.matrix(calc.statOutcomeByCovariate(field, value, quantile, outcome))
  colnames(quants) <- c("pathID", paste(c("0", "25", "50", "75", "100"), value, sep="_"))
  x <- merge(x, quants, by="pathID", all=TRUE)
  names(x) <- sub(outcome, value, names(x))

  form <- as.formula(paste(outcome, "~ pathID"))
  n <- aggregate(form, data=data.joined[data.joined[,field]==value,], length)
  x <- merge(x, n, by="pathID", all=TRUE)
  names(x)[length(names(x))] <- paste("N", value, sep="_")

  out <- x[x$N>cutoff,]

  return(out)
}


#' join.stratifiedOutcomeTables
#'
#' \code{join.stratifiedOutcomeTables}
#'
#' @param data
#' @param outcome
#' @return res.all
#'

join.stratifiedOutcomeTables <- function(data, outcome){
  ## res.all <- join.stratifiedOutcomeTables(data, outcome="start.to.diag")
  ## res.all <- join.stratifiedOutcomeTables(data, outcome="totalcost")
  ## write.csv(res.all, file="../../../output_data/pathway_summaries.csv")

  require(plyr)

  pathdata <- create.indivpathmatrix(data)
  data <- data.frame(data, pathdata)
  pathdata <- get.uniqueDataTestCombinations(pathdata)   #summary(pathdata)
  data.joined <- merge(data, pathdata)
  data.joined$dummy <- "all"
  data.joined[,outcome][data.joined[,outcome]<0] <- NA

  res.all <- join_all(dfs = list(get.meansdOutcomeByCovariate("dummy", "all", 0, outcome),
                                 get.meansdOutcomeByCovariate("DosanjhGrouped", "1", 0, outcome),
                                 get.meansdOutcomeByCovariate("DosanjhGrouped", "2", 0, outcome)),
# get.meansdOutcomeByCovariate("DosanjhGrouped", "3", 0, outcome),
# get.meansdOutcomeByCovariate("DosanjhGrouped", "4", 0, outcome)),
#                                  get.meansdOutcomeByCovariate("Sex", "Male", 0, outcome),
#                                  get.meansdOutcomeByCovariate("Sex", "Female", 0, outcome),
#                                  get.meansdOutcomeByCovariate("HIVpos", 0, 0, outcome),
#                                  get.meansdOutcomeByCovariate("HIVpos", 1, 0, outcome),
#                                  get.meansdOutcomeByCovariate("Ethnclass", "White", 0, outcome)),
                      type = "full", by = "pathID")

  res.all <- merge(pathdata, res.all)

  ## separate pathway and stats tables
  # merge(pathdata, get.meansdStartToDiagByCovariate("dummy", "all"))
  # merge(pathdata, get.meansdStartToDiagByCovariate("DosanjhGrouped", "1"))
  # merge(pathdata, get.meansdStartToDiagByCovariate("DosanjhGrouped", "2"))
  # merge(pathdata, get.meansdStartToDiagByCovariate("DosanjhGrouped", "3"))
  # merge(pathdata, get.meansdStartToDiagByCovariate("DosanjhGrouped", "4"))
  # merge(pathdata, get.meansdStartToDiagByCovariate("Sex", "Male"))
  # merge(pathdata, get.meansdStartToDiagByCovariate("Sex", "Female"))
  # merge(pathdata, get.meansdStartToDiagByCovariate("HIVpos", 0))
  # merge(pathdata, get.meansdStartToDiagByCovariate("HIVpos", 1))
  # merge(pathdata, get.meansdStartToDiagByCovariate("Ethnclass", "White"))

  form <- as.formula(paste(outcome, " ~ DosanjhGrouped"))
  print(aggregate(form, data=data.joined, summary, na.rm=T))
  form <- as.formula(paste(outcome, " ~ TBconfirmed"))
  print(aggregate(form, data=data.joined, summary, na.rm=T))

  round(data.matrix(res.all),2)
}


get.pooledtimeandcost <- function(){

##TODO##
}


#' Calculate summary statistics to include in final tables
#'
#' \code{maketable.Dosanjh_TimeCost}
#'
#' @param data
#'
#' @return matrix
#'
#' @seealso
#'
calc.summarystats_Dosanjh_TimeCost <- function(data){

  tab.time <- as.matrix(aggregate(start.to.diag~DosanjhGrouped, data=data, summary, na.rm=T))
  # tab <- cbind(tab, sd=as.matrix(aggregate(start.to.diag~DosanjhGrouped, data=data, sd, na.rm=T))[,2])
  tab.time <- rbind(tab.time, c("",summary(data$start.to.diag, na.rm=T)[-7])) #, sd=sd(data$start.to.diag, na.rm=T)))

  tab.cost <- as.matrix(aggregate(totalcost~DosanjhGrouped, data=data, summary, na.rm=T))
  # tab <- cbind(tab, sd=as.matrix(aggregate(start.to.diag~DosanjhGrouped, data=data, sd, na.rm=T))[,2])
  tab.cost <- rbind(tab.cost, c("",summary(data$totalcost, na.rm=T)[-7])) #, sd=sd(data$start.to.diag, na.rm=T)))

  return(cbind(tab.time[tab.time[,"DosanjhGrouped"]%in%c("1","2","3","4"),],
               tab.cost[tab.cost[,"DosanjhGrouped"]%in%c("1","2","3","4"),-1]))
}


#' Make test frequency table
#'
#' \code{make.testFreqTable}
#'
#' @param data
#' @return matrix

make.testFreqTable <- function(data){

  test.names <- c("TBcult", "QFN", "TSPOT", "TSTcut", "Smear", "BAL", #"CSF",
                  "HistBiop", "NeedleAsp", "PCR", "CXR", "CT", "MRI", "PET")

  data[,test.names] <- (data[,test.names]!="Not taken")

  nums <- aggregate(data[test.names], list(data$DosanjhGrouped), sum)
  names(nums) <- c("Dosanjh category", test.names)

  props <- as.matrix(round(nums[,test.names]/aggregate(data[test.names], list(data$DosanjhGrouped), length)[,test.names],2))
  nums <- as.matrix(nums)
  out <- matrix(paste(nums[,test.names], " (", props[,test.names], ")", sep=""), nrow=4)
  colnames(out) <- test.names
  out <- cbind(Dosanjh=c(1,2,3,4), out)
  out
}


#' Make test performance metrics table
#'
#' Includes true and false negatives and positive and negative predictive values.
#'
#' @param prop_highriskVECTOR Clinical judgement threshold as high risk patient
#' @param specificityVECTOR Rule-out test specificity
#' @return matrix

make.testPerformanceMetricsTable <- function(data,
                                             prop_highriskVECTOR = c(0.2, 0.4, 0.6),
                                             specificityVECTOR = c(0.9, 0.95, 0.99),
                                             sens = 0.9){
  out <- list()
  TN <- FN <- PPV <- NPV <- NNT <- NNH <- accuracy <- NULL

  for(i in 1:length(prop_highriskVECTOR)){

    out[[i]] <- make.ruleoutTable.pre(data=data, prop_highrisk = prop_highriskVECTOR[i])$combinedDosanjh

    for (j in 1:length(specificityVECTOR)){

      row <- out[[i]]$sensitivity==sens & out[[i]]$specificity==specificityVECTOR[j]

      TN <- c(TN, unique(out[[i]]$numTN[row]))
      FN <- c(FN, unique(out[[i]]$numFN[row]))
      PPV <- c(PPV, unique(out[[i]]$PPV[row]))
      NPV <- c(NPV, unique(out[[i]]$NPV[row]))
      NNT <- c(NNT, unique(out[[i]]$NNT[row]))
      NNH <- c(NNH, unique(out[[i]]$NNH[row]))
      accuracy <- c(accuracy, unique(out[[i]]$Accuracy[row]))
    }
  }

  out <- round(rbind(TN, FN, accuracy, PPV, NPV, NNT, NNH), 2)
  colnames(out) <- as.vector(outer(specificityVECTOR, prop_highriskVECTOR, paste))
  out
}


#' Make array of statistics for each bootstrap sample of IDEA study data
#'
#' @param n Number of bootstrap samples
#'
#' @return dim 3 array

boot.testPerformanceMetricTable <- function(n = 1000){

  ##TODO##
  #repeated this type of function for other functions (DRY)
  #consolidate inot a generic and curry

  dim0 <- dim(make.testPerformanceMetricsTable(data=data))
  barray <- array(NA, c(dim0, n))

  for(i in 1:n){
    bdata <- bootstrap.data(data)
    barray[,,i] <- make.testPerformanceMetricsTable(data=bdata)
  }

  meanvals <- apply(barray, c(1,2), function(x) round(mean(as.numeric(x)),2))
  sdvals <- apply(barray, c(1,2), function(x) round(sd(as.numeric(x)),2))

  return(matrix(paste(meanvals, " (",sdvals,")", sep=""), nrow = dim0[1]))
}


