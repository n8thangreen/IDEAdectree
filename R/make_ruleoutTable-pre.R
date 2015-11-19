#' make.ruleoutTable.pre
#'
#' \code{make.ruleoutTable.pre} creates the times and costs of a diagnositic pathway for
#' suspected active TB with and without an initial rule-out test, split by Dosanjh category.
#'
#' @param thresh test sensitivity and specificity
#' @param Ctest cost of rule-out test
#' @param FNcost false negative cost of true diagnosis or cost to start of standard pathway
#' @param FNtime false negative time to true diagnosis or start of pathway
#' @param ruleouttime time taken for rule-ou test result
#' @param pathreturn does a ruled-out individual return to the pathway? false positives return to standard pathway Y=1/N=0
#' @param qaly QALY for disease
#' @param A cost of one day in full health
#' @param npatients number of patients in cohort
#' @param cat4percent percent of patients in Dosanjh category 4
#' @param comb combined sensitivity, specificity and rule-out test cost array
#' @param cat3TB are Dosanjh category 3 patients all active TB?
#' @param cat4propfollowup proportion of negative Dosanjh category 4 patients, not immediately on standard pathway, who are followed-up at 6 weeks (alpha). We assume that all active TB cases are certain to be followed-up.
#' @param prop_highrisk Minimum predictive probability of TB risk score (delta)
#' @param               [previously, proportion of patients put on standard pathway by clinical judgment (gamma)],
#' @param stat Which statistic to use for time and cost estimate (e.g. Mean, Median, 1st Qu.)
#' @param model Which model structure/tree design to use (
#' \code{pretest.fixed} is include highest risk proportion with clinical judgement before rule-out test (fixed proportion risk factor independent) USED IN PAPER,
#' \code{pretest.var} is include highest risk proportion with clinical judgement before rule-out test (risk factor dependent),
#' \code{posttest.fixed} is test everyone first and include randomly selected proportion then include highest risk proportions,
#' \code{posttest.var} is test everyone first then remove highest risk proportion,
#' \code{pretest.var.sensspec.var} is remove highest risk proportion before rule-out test and modify sensitivity and specificity proportions wrt subset case-mix.)
#' @return list


make.ruleoutTable.pre <- function(
  # data = data,
  thresh = seq(from=1, to=0.8, by=-0.01),
  Ctest = c(200, 300),#, 200, 300, 400, 500, 600),
  FNcost = 0,
  FNtime = 42,  #6 weeks
  ruleouttime = 1L,
  pathreturn  = 1L,
  qaly = 0.67,
  A = 55,
  npatients = nrow(data),
  cat4percent = table(data$DosanjhGrouped)[4]*100/npatients,  #i.e no change
  comb=NA,
  cat3TB=TRUE,
  cat4propfollowup=0,
  prop_highrisk=0.4, #most of the fits are <0.7
  stat="Median",
  model="pretest.fixed"){


  ##TODO##
  #update prop_highrisk
  ## non-bootstrap sampled patients
  calc.costeqn.posttest.fixed <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    cost.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      (1-(1-prop_highrisk)*(1-cat4propfollowup)) * (pathreturn*pwaycost.old[cat]+FNcost)*numRuledOut.new[[cat]]
    }else (pathreturn*pwaycost.old[cat]+FNcost)*numRuledOut.new[[cat]]
    (pwaycost.old[cat]*numNotRuledOut.new[[cat]] + testcost*NumDosanjh[cat] + cost.FN)/NumDosanjh[cat]
  }

  calc.costeqn.posttest.var <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    cost.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      (1-(1-unlist(prop_highrisk))*(1-cat4propfollowup)) * (pathreturn*pwaycost.old[cat]+FNcost)*numRuledOut.new[[cat]]
    }else (pathreturn*pwaycost.old[cat]+FNcost)*numRuledOut.new[[cat]]
    (pwaycost.old[cat]*numNotRuledOut.new[[cat]] + testcost*NumDosanjh[cat] + cost.FN)/NumDosanjh[cat]
  }

  calc.costeqn.pretest.fixed <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    cost.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      pwaycost.old[cat]*numNotRuledOut.new[[cat]]
    }else (pwaycost.old[cat]*NumDosanjh[cat])
    ((prop_highrisk*pwaycost.old[cat]*NumDosanjh[cat]) + (1-prop_highrisk)*(testcost*NumDosanjh[cat] + cost.FN))/NumDosanjh[cat]
  }

  calc.costeqn.pretest.var <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    cost.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      pwaycost.lowrisk[cat]*numNotRuledOut.new[[cat]]
    }else (pwaycost.lowrisk[cat]*NumDosanjh[cat])
    ((prop_highrisk*pwaycost.highrisk[cat]*NumDosanjh[cat]) + (1-prop_highrisk)*(testcost*NumDosanjh[cat] + cost.FN))/NumDosanjh[cat]
  }

  ##---

  calc.timeToDiageqn.posttest.fixed <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    time.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      numRuledOut.new[[cat]] * (((1-(1-prop_highrisk)*(1-cat4propfollowup))*(pathreturn*daysToDiag.old[cat])) + ((1-prop_highrisk)*cat4propfollowup*FNtime))
    }else numRuledOut.new[[cat]]*((pathreturn*daysToDiag.old[cat]) + ((1-prop_highrisk)*FNtime))
    (daysToDiag.old[cat]*numNotRuledOut.new[[cat]] + ruleouttime*NumDosanjh[cat] + time.FN)/NumDosanjh[cat]
  }

  calc.timeToDiageqn.posttest.var <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    time.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      numRuledOut.new[[cat]] * (((1-(1-unlist(prop_highrisk))*(1-cat4propfollowup))*(pathreturn*daysToDiag.old[cat])) +
                                  ((1-unlist(prop_highrisk))*cat4propfollowup*FNtime))
    }else numRuledOut.new[[cat]]*((pathreturn*daysToDiag.old[cat]) + ((1-unlist(prop_highrisk))*FNtime))
    (daysToDiag.old[cat]*numNotRuledOut.new[[cat]] + ruleouttime*NumDosanjh[cat] + time.FN)/NumDosanjh[cat]
  }

  calc.timeToDiageqn.pretest.fixed <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    time.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      daysToDiag.old[cat]*numNotRuledOut.new[[cat]]
    }else (daysToDiag.old[cat]*NumDosanjh[cat] + numRuledOut.new[[cat]]*FNtime)
    ((prop_highrisk*daysToDiag.old[cat]*NumDosanjh[cat]) + (1-prop_highrisk)*(ruleouttime*NumDosanjh[cat] + time.FN))/NumDosanjh[cat]
  }

  calc.timeToDiageqn.pretest.var <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    time.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      daysToDiag.lowrisk[cat]*numNotRuledOut.new[[cat]]
    }else (daysToDiag.lowrisk[cat]*NumDosanjh[cat] + numRuledOut.new[[cat]]*FNtime)
    ((prop_highrisk*daysToDiag.highrisk[cat]*NumDosanjh[cat]) + (1-prop_highrisk)*(ruleouttime*NumDosanjh[cat] + time.FN))/NumDosanjh[cat]
  }


  calcCostEqn <- function(model, cat, cat3TB, cat4propfollowup, prop_highrisk){
    if(model=="pretest.fixed"){
      return(calc.costeqn.pretest.fixed(cat, cat3TB, cat4propfollowup, prop_highrisk))
    }else if(model=="pretest.var"){
      return(calc.costeqn.pretest.var(cat, cat3TB, cat4propfollowup, prop_highrisk))
    }else if(model=="posttest.fixed"){
      return(calc.costeqn.posttest.fixed(cat, cat3TB, cat4propfollowup, prop_highrisk))
    }else if(model=="posttest.var"){
      return(calc.costeqn.posttest.var(cat, cat3TB, cat4propfollowup, prop_highrisk))
    }else{stop("undefined model")}
  }

  calcTimeToDiagEqn <- function(model, cat, cat3TB, cat4propfollowup, prop_highrisk){
    if(model=="pretest.fixed"){
      return(calc.timeToDiageqn.pretest.fixed(cat, cat3TB, cat4propfollowup, prop_highrisk))
    }else if(model=="pretest.var"){
      return(calc.timeToDiageqn.pretest.var(cat, cat3TB, cat4propfollowup, prop_highrisk))
    }else if(model=="posttest.fixed"){
      return(calc.timeToDiageqn.posttest.fixed(cat, cat3TB, cat4propfollowup, prop_highrisk))
    }else if(model=="posttest.var"){
      return(calc.timeToDiageqn.posttest.var(cat, cat3TB, cat4propfollowup, prop_highrisk))
    }else{stop("undefined model")}
  }


  out <- list()
  nthresh <- length(thresh)

  NumDosanjh <- table(data$DosanjhGrouped)

  ## redistribute to change prevalence
  NumDosanjh[4] <- npatients/100*cat4percent
  NumDosanjh[1:3] <- (npatients-NumDosanjh[4])*prop.table(NumDosanjh[1:3])  #at the moment assumes that cat 3 are TB. dont think it makes much difference


  if(is.na(comb)){
    comb <- expand.grid(spec=thresh, sens=thresh, testcost=Ctest)
  }
  specificity <- comb$spec  #rep(thresh, each=length(nthresh))
  sensitivity <- comb$sens  #rep(thresh, by=nthresh)
  testcost <- comb$testcost

  ECDF.TB <- ecdf(data$riskfacScore[data$DosanjhGrouped%in%c(1,2)])
  ECDF.nonTB <- ecdf(data$riskfacScore[data$DosanjhGrouped==4])

  ## equivalent clinical judgement as a test dependent on risk factors
  spec.clinical <- ECDF.nonTB(prop_highrisk)
  sens.clinical <- 1-ECDF.TB(prop_highrisk)

  ## combined clinical judgement and rule-out test as new single test performance
  if(model=="posttest.var"){
    specificity <- pmin(specificity, spec.clinical)
    sensitivity <- pmax(sensitivity, sens.clinical)
    prop_highrisk.orig <- prop_highrisk
    prop_highrisk <- 0  #ie /delta=1
    spec.clinical <- 1
    sens.clinical <- 0
    model <- "posttest.fixed"
  }else if(model=="pretest.var.sensspec.var"){
    sensitivity <- pmax(0, 1 - ((1-sensitivity)/(1-sens.clinical)))
    specificity <- pmin(1, specificity/spec.clinical)
    model <- "pretest.var"
  }


  numRuledOut.new <- list()
  numRuledOut.new[[1]] <- NumDosanjh[1]*(1-sensitivity)
  numRuledOut.new[[2]] <- NumDosanjh[2]*(1-sensitivity)
  numRuledOut.new[[4]] <- NumDosanjh[4]*specificity

  if(cat3TB){numRuledOut.new[[3]] <- NumDosanjh[3]*(1-sensitivity)  #category 3 as active TB
  }else{numRuledOut.new[[3]] <- NumDosanjh[3]*specificity}          #category 3 as not active TB


  if(model=="posttest.var"){
#
#     ##TODO##
#     ## doesnt quite workout. shouldnt have a second kink in the INMB curve...
#
#     prop_highriskDosanjh <- list()
#
#     data.ordered <- data[order(data$riskfacScore),]
#
#     data1 <- data.ordered[data.ordered$DosanjhGrouped==1,]
#     data.post <- sapply(sensitivity, function(x) data1[1:(nrow(data1)*(1-x)),"riskfacScore"])   ##TODO## when x=1?
#     # data.post <- sapply(sensitivity, function(x) ifelse(x==1,NA,data1[1:(nrow(data1)*(1-x)),"riskfacScore"]))
#     # data.post <- sapply(sensitivity, function(x) data1[seq(length.out=nrow(data1)*(1-x)), "riskfacScore"])
#     prop_highriskDosanjh[[1]] <- sapply(data.post, function(x) sum(x>=prop_highrisk, na.rm=T)/sum(!is.na(x)))
#
#     data2 <- data.ordered[data.ordered$DosanjhGrouped==2,]
#     data.post <- sapply(sensitivity, function(x) data2[1:(nrow(data2)*(1-x)),"riskfacScore"])
#     # data.post <- sapply(sensitivity, function(x) ifelse(x==1,NA,data2[1:(nrow(data2)*(1-x)),"riskfacScore"]))
#     # data.post <- sapply(sensitivity, function(x) data2[seq(length.out=nrow(data2)*(1-x)), "riskfacScore"])
#     prop_highriskDosanjh[[2]] <- sapply(data.post, function(x) sum(x>=prop_highrisk, na.rm=T)/sum(!is.na(x)))
#
#     data3 <- data.ordered[data.ordered$DosanjhGrouped==3,]
#     if(cat3TB){
#       data.post <- sapply(sensitivity, function(x) data3[1:(nrow(data3)*(1-x)),"riskfacScore"])
#       # data.post <- sapply(sensitivity, function(x) ifelse(x==1,NA,data3[1:(nrow(data3)*(1-x)),"riskfacScore"]))
#       # data.post <- sapply(sensitivity, function(x) data3[seq(length.out=nrow(data3)*(1-x)),"riskfacScore"])
#     }else{
#       data.post <- sapply(specificity, function(x) data3[1:(nrow(data3)*x),"riskfacScore"])}
#     # data.post <- sapply(specificity, function(x) ifelse(x==0,NA,data3[1:(nrow(data3)*x),"riskfacScore"]))}
#       # data.post <- sapply(specificity, function(x) data3[seq(length.out=nrow(data3)*x), "riskfacScore"])}
#     prop_highriskDosanjh[[3]] <- sapply(data.post, function(x) sum(x>=prop_highrisk, na.rm=T)/sum(!is.na(x)))
#
#     data4 <- data.ordered[data.ordered$DosanjhGrouped==4,]
#     data.post <- sapply(specificity, function(x) data4[1:(nrow(data4)*x),"riskfacScore"])
#     # data.post <- sapply(specificity, function(x) ifelse(x==0,NA,data4[1:(nrow(data4)*x),"riskfacScore"]))
#     # data.post <- sapply(specificity, function(x) data4[seq(length.out=nrow(data4)*x), "riskfacScore"])
#     prop_highriskDosanjh[[4]] <- sapply(data.post, function(x) sum(x>=prop_highrisk, na.rm=T)/sum(!is.na(x)))

  }else{
    highriskPatients <- data$riskfacScore>=prop_highrisk

    prop_highriskDosanjh <- sum(data$DosanjhGrouped==1 & highriskPatients, na.rm=T)/sum(data$DosanjhGrouped==1 & !is.na(highriskPatients), na.rm=T)
    prop_highriskDosanjh[2] <- sum(data$DosanjhGrouped==2 & highriskPatients, na.rm=T)/sum(data$DosanjhGrouped==2 & !is.na(highriskPatients), na.rm=T)
    prop_highriskDosanjh[3] <- sum(data$DosanjhGrouped==3 & highriskPatients, na.rm=T)/sum(data$DosanjhGrouped==3 & !is.na(highriskPatients), na.rm=T)
    prop_highriskDosanjh[4] <- sum(data$DosanjhGrouped==4 & highriskPatients, na.rm=T)/sum(data$DosanjhGrouped==4 & !is.na(highriskPatients), na.rm=T)
  }


  numRuledOut.new.total <- rowSums(data.frame(numRuledOut.new[[1]],
                                              numRuledOut.new[[2]],
                                              numRuledOut.new[[3]],
                                              numRuledOut.new[[4]]), na.rm=TRUE)

  numNotRuledOut.new <- list()
  numNotRuledOut.new[[1]] <- NumDosanjh[1] - numRuledOut.new[[1]]
  numNotRuledOut.new[[2]] <- NumDosanjh[2] - numRuledOut.new[[2]]
  numNotRuledOut.new[[3]] <- NumDosanjh[3] - numRuledOut.new[[3]]
  numNotRuledOut.new[[4]] <- NumDosanjh[4] - numRuledOut.new[[4]]


  pwaycost.old    <- summary(data$totalcost[data$DosanjhGrouped=="1"])[stat]
  pwaycost.old[2] <- summary(data$totalcost[data$DosanjhGrouped=="2"])[stat]
  pwaycost.old[3] <- summary(data$totalcost[data$DosanjhGrouped=="3"])[stat]
  pwaycost.old[4] <- summary(data$totalcost[data$DosanjhGrouped=="4"])[stat]

  pwaycost.highrisk    <- summary(data$totalcost[data$DosanjhGrouped=="1" & highriskPatients])[stat]
  pwaycost.highrisk[2] <- summary(data$totalcost[data$DosanjhGrouped=="2" & highriskPatients])[stat]
  pwaycost.highrisk[3] <- summary(data$totalcost[data$DosanjhGrouped=="3" & highriskPatients])[stat]
  pwaycost.highrisk[4] <- summary(data$totalcost[data$DosanjhGrouped=="4" & highriskPatients])[stat]

  pwaycost.lowrisk    <- summary(data$totalcost[data$DosanjhGrouped=="1" & !highriskPatients])[stat]
  pwaycost.lowrisk[2] <- summary(data$totalcost[data$DosanjhGrouped=="2" & !highriskPatients])[stat]
  pwaycost.lowrisk[3] <- summary(data$totalcost[data$DosanjhGrouped=="3" & !highriskPatients])[stat]
  pwaycost.lowrisk[4] <- summary(data$totalcost[data$DosanjhGrouped=="4" & !highriskPatients])[stat]


  ##################
  # financial cost #
  ##################

  cost.new.notsampled <- list()

  # cost.new.notsampled <- plyr::alply(1:4, 1, function(x) calcCostEqn(model, x, cat3TB, cat4propfollowup, prop_highrisk))#, envir=parent.frame)  #fixed gamma
  cost.new.notsampled <- plyr::alply(1:4, 1, function(x) calcCostEqn(model, x, cat3TB, cat4propfollowup, prop_highriskDosanjh[x]))

  cost.new.notsampled.total <- rowSums(data.frame(cost.new.notsampled[[1]]*NumDosanjh[1],
                                                  cost.new.notsampled[[2]]*NumDosanjh[2],
                                                  cost.new.notsampled[[3]]*NumDosanjh[3],
                                                  cost.new.notsampled[[4]]*NumDosanjh[4]), na.rm=TRUE)/npatients

  diff_cost.old.new.notsampled <- list()
  diff_cost.old.new.notsampled[[1]]  <- pwaycost.old[1] - cost.new.notsampled[[1]]
  diff_cost.old.new.notsampled[[2]]  <- pwaycost.old[2] - cost.new.notsampled[[2]]
  diff_cost.old.new.notsampled[[3]]  <- pwaycost.old[3] - cost.new.notsampled[[3]]
  diff_cost.old.new.notsampled[[4]]  <- pwaycost.old[4] - cost.new.notsampled[[4]]
  diff_cost.old.new.notsampled.total <- rowSums(data.frame(diff_cost.old.new.notsampled[[1]]*NumDosanjh[1],
                                                           diff_cost.old.new.notsampled[[2]]*NumDosanjh[2],
                                                           diff_cost.old.new.notsampled[[3]]*NumDosanjh[3],
                                                           diff_cost.old.new.notsampled[[4]]*NumDosanjh[4]), na.rm=TRUE)/npatients

  ########
  # time #
  ########

  daysToDiag.old    <- summary(data$start.to.diag[data$DosanjhGrouped=="1"], na.rm=T)[stat]
  daysToDiag.old[2] <- summary(data$start.to.diag[data$DosanjhGrouped=="2"], na.rm=T)[stat]
  daysToDiag.old[3] <- summary(data$start.to.diag[data$DosanjhGrouped=="3"], na.rm=T)[stat]
  daysToDiag.old[4] <- summary(data$start.to.diag[data$DosanjhGrouped=="4"], na.rm=T)[stat]

  daysToDiag.highrisk    <- summary(data$start.to.diag[data$DosanjhGrouped=="1" & highriskPatients], na.rm=T)[stat]
  daysToDiag.highrisk[2] <- summary(data$start.to.diag[data$DosanjhGrouped=="2" & highriskPatients], na.rm=T)[stat]
  daysToDiag.highrisk[3] <- summary(data$start.to.diag[data$DosanjhGrouped=="3" & highriskPatients], na.rm=T)[stat]
  daysToDiag.highrisk[4] <- summary(data$start.to.diag[data$DosanjhGrouped=="4" & highriskPatients], na.rm=T)[stat]

  daysToDiag.lowrisk    <- summary(data$start.to.diag[data$DosanjhGrouped=="1" & !highriskPatients], na.rm=T)[stat]
  daysToDiag.lowrisk[2] <- summary(data$start.to.diag[data$DosanjhGrouped=="2" & !highriskPatients], na.rm=T)[stat]
  daysToDiag.lowrisk[3] <- summary(data$start.to.diag[data$DosanjhGrouped=="3" & !highriskPatients], na.rm=T)[stat]
  daysToDiag.lowrisk[4] <- summary(data$start.to.diag[data$DosanjhGrouped=="4" & !highriskPatients], na.rm=T)[stat]


  daysToDiag.new.notsampled <- list()

  # daysToDiag.new.notsampled <- plyr::alply(1:4, 1, function(x) calcTimeToDiagEqn(model,x,cat3TB,cat4propfollowup,prop_highrisk))  #fixed gamma
  daysToDiag.new.notsampled <- plyr::alply(1:4, 1, function(x) calcTimeToDiagEqn(model, x, cat3TB, cat4propfollowup, prop_highriskDosanjh[x]))
  daysToDiag.new.notsampled.total <- rowSums(data.frame(daysToDiag.new.notsampled[[1]]*NumDosanjh[1],
                                                        daysToDiag.new.notsampled[[2]]*NumDosanjh[2],
                                                        daysToDiag.new.notsampled[[3]]*NumDosanjh[3],
                                                        daysToDiag.new.notsampled[[4]]*NumDosanjh[4]), na.rm=TRUE)/npatients

  diff_daysToDiag.notsampled <- list()
  diff_daysToDiag.notsampled[[1]] <- daysToDiag.old[1] - daysToDiag.new.notsampled[[1]]
  diff_daysToDiag.notsampled[[2]] <- daysToDiag.old[2] - daysToDiag.new.notsampled[[2]]
  diff_daysToDiag.notsampled[[3]] <- daysToDiag.old[3] - daysToDiag.new.notsampled[[3]]
  diff_daysToDiag.notsampled[[4]] <- daysToDiag.old[4] - daysToDiag.new.notsampled[[4]]
  diff_daysToDiag.notsampled.total <- rowSums(data.frame(diff_daysToDiag.notsampled[[1]]*NumDosanjh[1],
                                                         diff_daysToDiag.notsampled[[2]]*NumDosanjh[2],
                                                         diff_daysToDiag.notsampled[[3]]*NumDosanjh[3],
                                                         diff_daysToDiag.notsampled[[4]]*NumDosanjh[4]), na.rm=TRUE)/npatients

  QALY.diff_time_total <- qaly*A*diff_daysToDiag.notsampled.total

  diff_total <- diff_cost.old.new.notsampled.total + QALY.diff_time_total


  ## more aggregated (i.e. not recorded by Dosanjh) than above ##

  if(model=="pretest.fixed"){

    Nnew <<- c(1-sens.clinical, 1-sens.clinical, ifelse(cat3TB,1-sens.clinical,spec.clinical), spec.clinical)%*%NumDosanjh
    Nnew_TB  <-  c(1-sens.clinical, 1-sens.clinical, ifelse(cat3TB,1-sens.clinical,0), 0)%*%NumDosanjh
    Nnew_nonTB <- c(0, 0, ifelse(cat3TB,0,spec.clinical), spec.clinical)%*%NumDosanjh
    prevNew <- Nnew_TB/Nnew


    testcostavoid <- spec.clinical*numRuledOut.new[[4]]*pwaycost.old[4]
    + ifelse(cat3TB, 0, spec.clinical*numRuledOut.new[[3]]*pwaycost.old[3])

    timeavoid <- spec.clinical*(numRuledOut.new[[4]]*daysToDiag.old[4] - NumDosanjh[4]*ruleouttime)
    + ifelse(cat3TB, 0, spec.clinical*(numRuledOut.new[[3]]*daysToDiag.old[3] - NumDosanjh[3]*ruleouttime))

    # testcostincur <- (1-prop_highriskDosanjh)%*%NumDosanjh *testcost
    testcostincur <- Nnew*testcost

    timeincur <- ((1-sens.clinical)*numRuledOut.new[[1]] + (1-sens.clinical)*numRuledOut.new[[2]])*FNtime +
      ((1-sens.clinical)*NumDosanjh[1] + (1-sens.clinical)*NumDosanjh[2])*ruleouttime +
      ifelse(cat3TB, (1-sens.clinical)*(numRuledOut.new[[3]]*FNtime + NumDosanjh[3]*ruleouttime), 0)

    calcCruleout_hat <- function(totalcostavoid, timecostincur, npatients)  (totalcostavoid - timecostincur)/Nnew


    ##TODO##
    riskprofile <- NA

    ## performance measures --
    ##https://en.wikipedia.org/wiki/Positive_and_negative_predictive_values

    ## positive/negative predictive value
    PPV <- sensitivity*prevNew/(sensitivity*prevNew + (1-specificity)*(1-prevNew))
    NPV <- specificity*(1-prevNew)/((1-sensitivity)*prevNew + specificity*(1-prevNew))

    numTN <- specificity*Nnew_nonTB
    numFN <- (1-sensitivity)*Nnew_TB
    numTP <- Nnew_TB - numFN
    accuracy <- (numTN + numTP)/Nnew


  }else if(model=="pretest.var"){
    testcostavoid <- (1-prop_highriskDosanjh[4])*numRuledOut.new[[4]]*pwaycost.lowrisk[4]
    + ifelse(cat3TB, 0, (1-prop_highriskDosanjh[3])*numRuledOut.new[[3]]*pwaycost.lowrisk[3])

    timeavoid <- (1-prop_highriskDosanjh[4])*(numRuledOut.new[[4]]*daysToDiag.lowrisk[4] - NumDosanjh[4]*ruleouttime)
    + ifelse(cat3TB, 0, (1-prop_highriskDosanjh[3])*(numRuledOut.new[[3]]*daysToDiag.lowrisk[3] - NumDosanjh[3]*ruleouttime))

    testcostincur <- ((1-prop_highriskDosanjh)%*%NumDosanjh) *testcost

    timeincur <- ((1-prop_highriskDosanjh[1])*numRuledOut.new[[1]] + (1-prop_highriskDosanjh[2])*numRuledOut.new[[2]])*FNtime +
      ((1-prop_highriskDosanjh[1])*NumDosanjh[1] + (1-prop_highriskDosanjh[2])*NumDosanjh[2])*ruleouttime +
      ifelse(cat3TB, (1-prop_highriskDosanjh[3])*(numRuledOut.new[[3]]*FNtime + NumDosanjh[3]*ruleouttime), 0)

    calcCruleout_hat <- function(totalcostavoid, timecostincur, npatients) (totalcostavoid-timecostincur)/((1-prop_highriskDosanjh)%*%NumDosanjh)


    riskprofile <- list(costs=data.frame(testcost + A*qaly*(ruleouttime+FNtime),
                                         testcost + A*qaly*ruleouttime,
                                         testcost + A*qaly*ruleouttime - pwaycost.old[4],
                                         0),
                        probs=data.frame((1-prop_highriskDosanjh[1])*numRuledOut.new[[1]]/NumDosanjh[1],
                                         (1-prop_highriskDosanjh[1])*numNotRuledOut.new[[1]]/NumDosanjh[1] + (1-prop_highriskDosanjh[4])*numNotRuledOut.new[[4]]/NumDosanjh[4],
                                         (1-prop_highriskDosanjh[4])*numRuledOut.new[[4]]/NumDosanjh[4],
                                         (NumDosanjh[1]*prop_highriskDosanjh[1] + NumDosanjh[4]*prop_highriskDosanjh[4])/(NumDosanjh[1]+NumDosanjh[4]) ))

  }else if(model=="posttest.fixed"){

    Nnew <<- npatients

    testcostavoid <- spec.clinical*(1-cat4propfollowup)*
      (numRuledOut.new[[4]]*pwaycost.old[4] + ifelse(cat3TB,0, numRuledOut.new[[3]]*pwaycost.old[3]))

    timeavoid <- numRuledOut.new[[4]]*spec.clinical*((1-cat4propfollowup)*daysToDiag.old[4] - (cat4propfollowup*FNtime)) -
      (NumDosanjh[4]*ruleouttime) +
      ifelse(cat3TB,0,
             numRuledOut.new[[3]]*spec.clinical*((1-cat4propfollowup)*daysToDiag.old[3] - (cat4propfollowup*FNtime)) - (NumDosanjh[3]*ruleouttime))

    testcostincur <- testcost*npatients

    timeincur <- (numRuledOut.new[[1]]+numRuledOut.new[[2]])*FNtime*(1-sens.clinical) + (sum(NumDosanjh[1:2])*ruleouttime) +
      ifelse(cat3TB, numRuledOut.new[[3]]*FNtime*(1-sens.clinical) + (NumDosanjh[3]*ruleouttime), 0)


    calcCruleout_hat <- function(totalcostavoid, timecostincur, npatients) (totalcostavoid-timecostincur)/npatients

    riskprofile <- NA   ##TODO##
    PPV <- NPV <- NA
    TN <- FN <- NA
    accuracy <- NA

  }else if(model=="posttest.var"){
    ##TODO##
    # use lowrisk times/costs?

    testcostavoid <- (1-prop_highriskDosanjh[[4]])*(1-cat4propfollowup)*numRuledOut.new[[4]]*pwaycost.old[4] +
      ifelse(cat3TB,0, (1-prop_highriskDosanjh[[3]])*(1-cat4propfollowup)*numRuledOut.new[[3]]*pwaycost.old[3])

    timeavoid <- numRuledOut.new[[4]]*(1-prop_highriskDosanjh[[4]])*((1-cat4propfollowup)*daysToDiag.old[4] - (cat4propfollowup*FNtime)) -
      (NumDosanjh[4]*ruleouttime) +
      ifelse(cat3TB,0,
             numRuledOut.new[[3]]*(1-prop_highriskDosanjh[[3]])*((1-cat4propfollowup)*daysToDiag.old[3] - (cat4propfollowup*FNtime)) -
               (NumDosanjh[3]*ruleouttime))

    testcostincur <- testcost*npatients

    timeincur <- numRuledOut.new[[1]]*FNtime*(1-prop_highriskDosanjh[[1]]) +
      numRuledOut.new[[2]]*FNtime*(1-prop_highriskDosanjh[[2]]) +
      (sum(NumDosanjh[1:2])*ruleouttime) +
      ifelse(cat3TB, numRuledOut.new[[3]]*FNtime*(1-prop_highriskDosanjh[[3]]) + (NumDosanjh[3]*ruleouttime), 0)

    calcCruleout_hat <- function(totalcostavoid, timecostincur, npatients) (totalcostavoid-timecostincur)/npatients
  }




  timecostavoid  <- A*qaly*timeavoid
  timecostincur  <- A*qaly*timeincur
  totalcostavoid <- testcostavoid + timecostavoid
  totalcostincur <- testcostincur + timecostincur
  totalcostdiff  <- totalcostavoid - totalcostincur

  stdpway.tcost.Dsnjh.LOW  <- pwaycost.lowrisk + A*qaly*daysToDiag.lowrisk
  stdpway.tcost.Dsnjh.HIGH <- pwaycost.highrisk + A*qaly*daysToDiag.highrisk
  mean.stdpway.tcost.HIGH.TB <- mean(stdpway.tcost.Dsnjh.HIGH[1:3])
  mean.stdpway.tcost.LOW.TB  <- mean(stdpway.tcost.Dsnjh.LOW[1:3])
  ruleouttest.totalcost <- testcost + A*qaly*ruleouttime

  Cruleout_hat <- calcCruleout_hat(totalcostavoid, timecostincur, npatients)
  INMB <- totalcostdiff/npatients
  INHB <- qaly*(timeavoid-timeincur)/(365*npatients) - (testcostincur-testcostavoid)/(365*A*npatients)
  ICER <- (testcostavoid-testcostincur)/(qaly*(timeavoid-timeincur))


  ## cost-neutral test cost summary stats
  ##TODO## pretest.fixed only
  sens1 <- spec1 <- 0.90
  C1 <- round(unique(Cruleout_hat[sensitivity==sens1 & specificity==spec1]),0)
  C_epsilon <- 20
  C2 <- C1 + C_epsilon

  C2.bestmatch.sens <- abs(Cruleout_hat[specificity==spec1]-C2) == min(abs(Cruleout_hat[specificity==spec1]-C2))
  C2.bestmatch.spec <- abs(Cruleout_hat[sensitivity==sens1]-C2) == min(abs(Cruleout_hat[sensitivity==sens1]-C2))

  sens2 <- unique(sensitivity[specificity==spec1][C2.bestmatch.sens])
  spec2 <- unique(specificity[sensitivity==sens1][C2.bestmatch.spec])

  C_hat.sensSpecRatio.exact <- (Nnew_TB*A*qaly*FNtime)/(Nnew_nonTB*(pwaycost.old[4]+A*qaly*daysToDiag.old[4]))   ##TODO## test
  C_hat.sensSpecRatio.numerical <- (spec2-spec1)/(sens2-sens1)


  deltaSpec <- Nnew/(Nnew_nonTB*(pwaycost.old[4]+A*qaly*daysToDiag.old[4]))
  deltaSens <- Nnew/(Nnew_TB*A*qaly*FNtime)
  deltaHypotenuse <- sqrt(deltaSpec^2 + deltaSens^2)

  C_hat.sensSpecSlope.exact <- deltaSpec*deltaSens/deltaHypotenuse
  C_hat.sensSpecSlope.numerical  <- ((spec2-spec1)*(sens2-sens1)/(sqrt((spec2-spec1)^2 + (sens2-sens1)^2)))/C_epsilon


  #because change raw values earlier to _effective_ values
  #for var models
  specificity.effective <- specificity
  sensitivity.effective <- sensitivity
  specificity <- comb$spec  #rep(thresh, each=length(nthresh))
  sensitivity <- comb$sens  #rep(thresh, by=nthresh)
  testcost <- comb$testcost


  out.combined <- data.frame(testcost=testcost,
                             sensitivity=sensitivity,
                             specificity=specificity,
                             sensitivity.effective=sensitivity.effective,
                             specificity.effective=specificity.effective,
                             testcostavoid,
                             timeavoid,
                             timecostavoid,
                             testcostincur,
                             timeincur,
                             timecostincur,
                             totalcostavoid,
                             totalcostincur,
                             totalcostdiff,
                             INMB,
                             INHB,
                             Cruleout_hat,
                             ICER,
                             PPV=PPV,
                             NPV=NPV,
                             numTN=numTN,
                             numFN=numFN,
                             Accuracy=accuracy)

  surfacevals <- c(C_hat.sensSpecRatio.exact = C_hat.sensSpecRatio.exact,
                   C_hat.sensSpecRatio.numerical = C_hat.sensSpecRatio.numerical,
                   C_hat.sensSpecSlope.exact = C_hat.sensSpecSlope.exact,
                   C_hat.sensSpecSlope.numerical = C_hat.sensSpecSlope.numerical)


  fixedcosts <- c(stdpway.cost.Dsnjh.LOW =stdpway.tcost.Dsnjh.LOW,
                  stdpway.cost.Dsnjh.HIGH=stdpway.tcost.Dsnjh.HIGH)

  ##TODO##
  ## branching points expected values
  ## EV3 and prop_highriskDosanjh for all positive \gamma
  EV <- c(EV1 <- mean.stdpway.tcost.LOW.TB + (1-sensitivity)*A*qaly*FNtime,
          EV2 <- (1-specificity)*stdpway.tcost.Dsnjh.LOW[4],
          EV3 <- sum(unlist(prop_highriskDosanjh[1])*mean.stdpway.tcost.HIGH.TB, (1-unlist(prop_highriskDosanjh[1]))*(ruleouttest.totalcost + ifelse(is.na(EV1),0,EV1)), na.rm=T),
          EV4 <- sum(unlist(prop_highriskDosanjh[4])*stdpway.tcost.Dsnjh.HIGH[4], (1-unlist(prop_highriskDosanjh[4]))*(ruleouttest.totalcost + ifelse(is.na(EV2),0,EV2)), na.rm=T),
          ((npatients-NumDosanjh[4])*ifelse(is.na(EV3),0,EV3) + NumDosanjh[4]*ifelse(is.na(EV4),0,EV4))/npatients)
  names(EV) <- c("EV1", "EV2", "EV3", "EV4", "EV5")

  out.costbyDosanjh <- round(data.frame(#scenarioNum=1:length(sensitivity),
    testcost=testcost,
    #FNcost=FNcost,
    sensitivity=sensitivity,
    specificity=specificity,
    numRuledOut.new1=round(numRuledOut.new[[1]]),
    numRuledOut.new2=round(numRuledOut.new[[2]]),
    numRuledOut.new3=round(numRuledOut.new[[3]]),
    numRuledOut.new4=round(numRuledOut.new[[4]]),
    numRuledOut.new=round(numRuledOut.new.total),
    #numRuledOut.final=NumDosanjh[4],
    pwaycost.old1=pwaycost.old[1],
    pwaycost.old2=pwaycost.old[2],
    pwaycost.old3=pwaycost.old[3],
    pwaycost.old4=pwaycost.old[4],
    cost.new1=cost.new.notsampled[[1]],
    cost.new2=cost.new.notsampled[[2]],
    cost.new3=cost.new.notsampled[[3]],
    cost.new4=cost.new.notsampled[[4]],
    cost.new.total=cost.new.notsampled.total,
    diff_cost1=diff_cost.old.new.notsampled[[1]],
    diff_cost2=diff_cost.old.new.notsampled[[2]],
    diff_cost3=diff_cost.old.new.notsampled[[3]],
    diff_cost4=diff_cost.old.new.notsampled[[4]],
    diff_cost.total=diff_cost.old.new.notsampled.total), 2)

  out.diagtimebyDosanjh <- round(data.frame(#scenarioNum=1:length(sensitivity),
    sensitivity=sensitivity,
    specificity=specificity,
    #FNtime=FNtime,
    #ruleouttime=ruleouttime,
    numRuledOut.new1=round(numRuledOut.new[[1]]),
    numRuledOut.new2=round(numRuledOut.new[[2]]),
    numRuledOut.new3=round(numRuledOut.new[[3]]),
    numRuledOut.new4=round(numRuledOut.new[[4]]),
    numRuledOut.new=round(numRuledOut.new.total),
    #numRuledOut.final=NumDosanjh[4],
    #daysToDiag.old1=daysToDiag.old[1],
    #daysToDiag.old2=daysToDiag.old[2],
    #daysToDiag.old3=daysToDiag.old[3],
    #daysToDiag.old4=daysToDiag.old[4],
    daysToDiag.new1=daysToDiag.new.notsampled[[1]],
    daysToDiag.new2=daysToDiag.new.notsampled[[2]],
    daysToDiag.new3=daysToDiag.new.notsampled[[3]],
    daysToDiag.new4=daysToDiag.new.notsampled[[4]],
    daysToDiag.new.total=daysToDiag.new.notsampled.total,
    diff_daysToDiag1=diff_daysToDiag.notsampled[[1]],
    diff_daysToDiag2=diff_daysToDiag.notsampled[[2]],
    diff_daysToDiag3=diff_daysToDiag.notsampled[[3]],
    diff_daysToDiag4=diff_daysToDiag.notsampled[[4]],
    diff_daysToDiag.total=diff_daysToDiag.notsampled.total), 2)

  # write.csv(out.costbyDosanjh, "../../../output_data/ruleout-costtable-notsampled.csv")
  # write.csv(out.diagtimebyDosanjh, "../../../output_data/ruleout-diagtimetable-notsampled.csv")

  list(combinedDosanjh=out.combined, costbyDosanjh=out.costbyDosanjh, diagtimebyDosanjh=out.diagtimebyDosanjh, EV=EV, fixedcosts=fixedcosts,
       prop_highriskDosanjh=prop_highriskDosanjh, ruleout.totalcost=ruleouttest.totalcost,
       pwaycost.highrisk=pwaycost.highrisk, pwaycost.lowrisk=pwaycost.lowrisk,
       riskprofile=riskprofile, surfacevals=surfacevals)
}


