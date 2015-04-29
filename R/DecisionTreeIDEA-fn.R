#' relabelDecTreeOutcomes
#'
#' \code{relabelDecTreeOutcomes} returns joined tables of inputs and predicted outcomes
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

relabelDecTreeOutcomes <- function(dt, join_names){
  # append filtering label to column names

  for (i in names(dt)){
    a <- !(names(dt[[as.character(i)]]) %in% join_names)
    names(dt[[i]])[a] <- paste(names(dt[[i]])[a], "_", i, sep="")
  }
  dt
}


#' addMissingEntries
#'
#' \code{addMissingEntries}
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

addMissingEntries <- function(names, dt){
  # fill-in empty columns with zeros

  targetnames <- c(names,"1_TRUE","2_TRUE","3_TRUE","4_TRUE","1_FALSE","2_FALSE","3_FALSE","4_FALSE")
  missingnames <- targetnames[!targetnames%in%names(dt)]
  missingmatrix <- matrix(0, nrow=nrow(dt), ncol=length(missingnames))
  colnames(missingmatrix) <- missingnames
  dt <- data.frame(dt, missingmatrix, check.names=F)
  dt <- dt[,targetnames]

  len <- length(names)
  ncols <- ncol(dt)
  dt[,(len+1):ncols][is.na(dt[,(len+1):ncols])] <- 0
  dt
}


#' calcSumRow
#'
#' \code{calcSumRow}
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

calcSumRow <- function(data, test, subsetvar, namesdt, val){
  # sum down rows for each column

  if(missing(val)){
    val <- unique(data[,test])
    txt <- ""
  }else{
    txt <- val}

  dat <- data[data[,test]%in%val,]

  c(txt,"","","",
    tabulate(dat[dat[,subsetvar],"DosanjhGrouped"], nbins=4),
    tabulate(dat[!dat[,subsetvar],"DosanjhGrouped"], nbins=4))
}

#' createDecisionTree
#'
#' \code{createDecisionTree}
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

createDecisionTree <- function(data, test.DT, diag_obs.DT, diag_est.DT){
  # calculate deterministic decision pathway
  # for different combinations of estimated and
  # observed steps on the decision tree

  ##TODO##
  # should these be $data as well??

  out <- list()

  out$t1 <- joinDataAndRules(data, test.DT)

  mat_treat1est <- out$t1$data
  out$t1$data$treat1.est <- out$t1$data$treat1
  out$t1$data$treat1 <- out$t1$data$step1Diag   #because join on treat1
  out$t1$data <- out$t1$data[,names(out$t1$data)!="step1Diag"]

  out$treat1obs_Sxobs <- joinDataAndRules(out$t1$data, diag_obs.DT)
  out$treat1obs_Sxest <- joinDataAndRules(out$t1$data, diag_est.DT)
  out$treat1est_Sxest <- joinDataAndRules(mat_treat1est, diag_est.DT)

  out$treat1obs_Sxest$dt <- out$treat1est_Sxest$dt <- diag_obs.DT

  out$treat1obs_Sxest$data$"Sx_resolved.obs" <- out$treat1obs_Sxest$data$"Sx_resolved"
  out$treat1obs_Sxest$data$"Sx_resolved" <- out$treat1obs_Sxest$data$"Sx_resolved.est"
  out$treat1est_Sxest$data$"Sx_resolved.obs" <- out$treat1est_Sxest$data$"Sx_resolved"
  out$treat1est_Sxest$data$"Sx_resolved" <- out$treat1est_Sxest$data$"Sx_resolved.est"

  out
}

#' joinDataAndRules
#'
#' \code{joinDataAndRules}
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

joinDataAndRules <- function(data, DTobj){
  x <- list()
  x$data <- merge(data, DTobj$all, by=names(DTobj$input), all.x=T)
  x$data$Freq <- 1
  x$dt <- DTobj
  x
}


#' freqTable.DecisionTree
#'
#' \code{freqTable.DecisionTree}
#'
#' @param decntree
#' @param subsetvar
#' @return one combined frequency table of each decision outcomes

freqTable.DecisionTree <- function(decntree, subsetvar){
  #
  # type: test (pre culture) or diag (post culture)
  # call: freqTable.DecisionTree(out$treat1obs_Sxest, preCult=TRUE, "step1Diag")

  names <- names(decntree$dt$input)
  test <- names[1]
  formula <- as.formula(paste(paste(names, collapse="+"),"~", "DosanjhGrouped |", subsetvar))

  dtNEG <- decntree$data[decntree$data[,test]=="NEGATIVE",]
  dtNEG <- reshape::cast(dtNEG, formula, fun=sum, value="Freq")
  dtNEG <- relabelDecTreeOutcomes(dt=dtNEG, join_names=names)
  dt <- plyr::join_all(dtNEG, by=names, type="full")
  dt <- plyr::ddply(dt, names)  #order rows
  dt <- addMissingEntries(names, dt)

  dt <- rbind(calcSumRow(data=decntree$data, test, subsetvar, namesdt=names(dtNEG), val="POSITIVE"),
              # calcSumRow(decntree, test, subsetvar, "Not taken"),
              dt,
              calcSumRow(decntree$data, test, subsetvar, names(dtNEG)))

  dt <- data.frame(dt,  Total=apply(dt[,names(dt)!=names], 1, function(x) sum(as.numeric(x))), check.names = FALSE)

  dt
}


#' writeAllFreqTables
#'
#' \code{writeAllFreqTables}
#'
#' @param out joined individual data and idealised pathway
#' @param dir output folder
#' @return list of frequency tables


writeAllFreqTables <- function(out, dir){

  nm <- names(out)[names(out)!="t1"]

  ## before culture ##
  treat1obs <- freqTable.DecisionTree(decntree=out$t1, subsetvar="treat1")
  treat1est <- freqTable.DecisionTree(decntree=out$t1, subsetvar="treat1.est")
  write.csv(treat1obs, file=paste(dir,"/treat1obs.csv", sep=""))
  write.csv(treat1est, file=paste(dir,"/treat1est.csv", sep=""))

  ## after culture ##
  for (i in nm){
    assign(i, freqTable.DecisionTree(decntree=out[[i]], subsetvar="treat2.ideal"))  #could just use a temp var instead?
    file <- paste('"', dir, '/', i, '.csv")', sep="")
    eval(parse(text = paste('write.csv(', i, ', file=', file, sep="")))
  }

  res <- list()
  res$treat1obs <- treat1obs
  res$treat1est <- treat1est

  invisible(list(res, mget(nm)))
}


#' create DT object
#'
#' \code{createDTobj}
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

createDTobj <- function(rules, input){
  res <- list()
  res$input <- rules[,input, drop=FALSE]
  res$out <- rules[, names(rules)[!names(rules)%in%input], drop=FALSE]
  res$all <- rules
  res
}


#' createCombinedGrid
#'
#' \code{createCombinedGrid}
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

createCombinedGrid <- function(dir = "C:/Users/Nathan/Dropbox/TB/IDEA/output_data/idealpathway_grids/"){
  ##
  ## createCombinedGrid()

  wd.old <- getwd()
  setwd(dir)

  rules.tests <- read.csv(file="tests.csv")
  rules.diag.est <- read.csv(file="diag_Sxest.csv")
  rules.diag.obs <- read.csv(file="diag_Sxobs.csv")

  test.DT <- createDTobj(rules.tests, input = c("Smear","TSTres","CXR","Risk_factors"))
  diag_est.DT <- createDTobj(rules.diag.est, input = c("TBcult","treat1","TBconfirmed","Alt_diag"))
  diag_obs.DT <- createDTobj(rules.diag.obs, input = c("TBcult","treat1","Alt_diag","Sx_resolved"))

  save(test.DT, diag_est.DT, diag_obs.DT, file="allgrids.RData")
  setwd <- wd.old
}


