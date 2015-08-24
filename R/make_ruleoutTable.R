#' make.ruleoutTable
#'
#' \code{make.ruleoutTable} creates the times and costs of a diagnositic pathway for
#' suspected active TB with and without an initial rule-out test, split by Dosanjh category.
#'
#' @param thresh test sensitivity and specificity
#' @param Ctest cost of rule-out test
#' @param FNcost false negative cost
#' @param FNtime false negative time
#' @param ruleouttime time for rule-ou test result
#' @param pathreturn does a ruled-out individual return to the pathway?
#' @param qaly QALY for disease
#' @param A cost of one day in full health
#' @param npatients number of patients in cohort
#' @param cat4percent percent of patients in Dosanjh category 4
#' @param comb combined sensitivity, specificity and rule-out test cost array
#' @param cat3TB are Dosanjh category 3 patients all active TB?
#' @param cat4propfollowup proportion of negative Dosanjh category 4 patients, not immediately on standard pathway, who are followed-up at 6 weeks (alpha)
#' @param prop_highrisk proportion of negative active TB patients put immediately on standard pathway (gamma)
#' @return list


make.ruleoutTable <- function(
  thresh = seq(from=1, to=0.7, by=-0.01),   #test sensitivity and specificity
  Ctest = c(50, 100),
  FNcost = 0,   #false negative cost of true diagnosis or cost to start of standard pathway
  FNtime = 42L,  #6 weeks #false negative time to true diagnosis or start of pathway
  ruleouttime = 1L,
  pathreturn  = 1L, #false positives return to standard pathway Y=1/N=0
  qaly = 0.67,
  A = 55,
  npatients = nrow(data),
  cat4percent = table(data$DosanjhGrouped)[4]*100/npatients,  #i.e no change
  comb=NA,
  cat3TB=TRUE,
  cat4propfollowup=0,
  prop_highrisk=0.72,
  stat="Mean",
  model="highriskAtStart"){


  ## non-bootstrap sampled patients
  calc.costeqn.testFirst <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    cost.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      (1-(1-prop_highrisk)*(1-cat4propfollowup)) * (pathreturn*pwaycost.old[cat]+FNcost)*numRuledOut.new[[cat]]
    }else (pathreturn*pwaycost.old[cat]+FNcost)*numRuledOut.new[[cat]]
    (pwaycost.old[cat]*numNotRuledOut.new[[cat]] + testcost*NumDosanjh[cat] + cost.FN)/NumDosanjh[cat]
  }

  calc.costeqn.highriskAtStart <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    cost.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      pwaycost.old[cat]*numNotRuledOut.new[[cat]]
    }else (pwaycost.old[cat]*NumDosanjh[cat])
    ((prop_highrisk*pwaycost.old[cat]*NumDosanjh[cat]) + (1-prop_highrisk)*(testcost*NumDosanjh[cat] + cost.FN))/NumDosanjh[cat]
  }

  calc.timeToDiageqn.testFirst <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    time.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      numRuledOut.new[[cat]] * (((1-(1-prop_highrisk)*(1-cat4propfollowup))*(pathreturn*daysToDiag.old[cat])) + ((1-prop_highrisk)*cat4propfollowup*FNtime))
    }else numRuledOut.new[[cat]]*((pathreturn*daysToDiag.old[cat]) + ((1-prop_highrisk)*FNtime))
    (daysToDiag.old[cat]*numNotRuledOut.new[[cat]] + ruleouttime*NumDosanjh[cat] + time.FN)/NumDosanjh[cat]
  }

  calc.timeToDiageqn.highriskAtStart <- function(cat, cat3TB, cat4propfollowup, prop_highrisk){
    time.FN <- if(cat==4 | (cat==3 & cat3TB==F)){
      daysToDiag.old[cat]*numNotRuledOut.new[[cat]]
    }else (daysToDiag.old[cat]*NumDosanjh[cat] + numRuledOut.new[[cat]]*FNtime)
    ((prop_highrisk*daysToDiag.old[cat]*NumDosanjh[cat]) + (1-prop_highrisk)*(ruleouttime*NumDosanjh[cat] + time.FN))/NumDosanjh[cat]
  }


  calcCostEqn <- function(model, cat, cat3TB, cat4propfollowup, prop_highrisk){
    if(model=="highriskAtStart"){
      return(calc.costeqn.highriskAtStart(cat, cat3TB, cat4propfollowup, prop_highrisk))
    }else{
      return(calc.costeqn.testFirst(cat, cat3TB, cat4propfollowup, prop_highrisk))}
  }
  calcTimeToDiagEqn <- function(model, cat, cat3TB, cat4propfollowup, prop_highrisk){
    if(model=="highriskAtStart"){
      return(calc.timeToDiageqn.highriskAtStart(cat, cat3TB, cat4propfollowup, prop_highrisk))
    }else{
      return(calc.timeToDiageqn.testFirst(cat, cat3TB, cat4propfollowup, prop_highrisk))}
  }


  out <- list()
  nthresh <- length(thresh)

  NumDosanjh <- table(data$DosanjhGrouped)
  NumDosanjh[4] <- npatients/100*cat4percent
  NumDosanjh[1:3] <- (npatients-NumDosanjh[4])*prop.table(NumDosanjh[1:3])  #at the moment assumes that cat 3 are TB. dont think it makes much difference

  if(is.na(comb)){
    comb <- expand.grid(spec=thresh, sens=thresh, testcost=Ctest)
    # comb <- expand.grid( c(1, 0.7), c(1, 0.9, 0.8, 0.7), Ctest)   #Excel values
  }
  specificity <- comb$spec  #rep(thresh, each=length(nthresh))
  sensitivity <- comb$sens  #rep(thresh, by=nthresh)
  testcost <- comb$testcost

  numRuledOut.new  <- list()
  numRuledOut.new[[1]] <- NumDosanjh[1]*(1-sensitivity)
  numRuledOut.new[[2]] <- NumDosanjh[2]*(1-sensitivity)
  numRuledOut.new[[4]] <- NumDosanjh[4]*specificity

  if(cat3TB){numRuledOut.new[[3]] <- NumDosanjh[3]*(1-sensitivity)  #category 3 as active TB
  }else{numRuledOut.new[[3]] <- NumDosanjh[3]*specificity}          #category 3 as not active TB

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


  ##################
  # financial cost #
  ##################

  cost.new.notsampled <- list()

  # cost.new.notsampled[[1]] <- (numNotRuledOut.new[[1]]*pwaycost.old[1] + NumDosanjh[1]*testcost + numRuledOut.new[[1]]*pwaycost.old[1])/NumDosanjh[1]
  # cost.new.notsampled[[2]] <- (numNotRuledOut.new[[2]]*pwaycost.old[2] + NumDosanjh[2]*testcost + numRuledOut.new[[2]]*pwaycost.old[2])/NumDosanjh[2]
  # cost.new.notsampled[[3]] <- (numNotRuledOut.new[[3]]*pwaycost.old[3] + NumDosanjh[3]*testcost + numRuledOut.new[[3]]*pwaycost.old[3])/NumDosanjh[3]
  # cost.new.notsampled[[4]] <- (numNotRuledOut.new[[4]]*pwaycost.old[4] + NumDosanjh[4]*testcost)/NumDosanjh[4]

  cost.new.notsampled <- plyr::alply(1:4, 1, function(x) calcCostEqn(model,x,cat3TB,cat4propfollowup,prop_highrisk))#, envir=parent.frame)
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

  pwaycost.old    <- summary(data$start.to.diag[data$DosanjhGrouped=="1"], na.rm=T)[stat]
  pwaycost.old[2] <- summary(data$start.to.diag[data$DosanjhGrouped=="2"], na.rm=T)[stat]
  pwaycost.old[3] <- summary(data$start.to.diag[data$DosanjhGrouped=="3"], na.rm=T)[stat]
  pwaycost.old[4] <- summary(data$start.to.diag[data$DosanjhGrouped=="4"], na.rm=T)[stat]


  daysToDiag.new.notsampled <- list()
  # daysToDiag.new.notsampled[[1]] <- (numNotRuledOut.new[[1]]*daysToDiag.old[1] + numRuledOut.new[[1]]*(FNtime+daysToDiag.old[1]) + NumDosanjh[1]*ruleouttime)/NumDosanjh[1]
  # daysToDiag.new.notsampled[[2]] <- (numNotRuledOut.new[[2]]*daysToDiag.old[2] + numRuledOut.new[[2]]*(FNtime+daysToDiag.old[2]) + NumDosanjh[2]*ruleouttime)/NumDosanjh[2]
  # daysToDiag.new.notsampled[[3]] <- (numNotRuledOut.new[[3]]*daysToDiag.old[3] + numRuledOut.new[[3]]*(FNtime+daysToDiag.old[3]) + NumDosanjh[3]*ruleouttime)/NumDosanjh[3]
  # daysToDiag.new.notsampled[[4]] <- (numNotRuledOut.new[[4]]*daysToDiag.old[4] + NumDosanjh[4]*ruleouttime)/NumDosanjh[4]

  daysToDiag.new.notsampled <- plyr::alply(1:4, 1, function(x) calcTimeToDiagEqn(model,x,cat3TB,cat4propfollowup,prop_highrisk))
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


  ## more aggregated (not recorded by Dosanjh) than above
  ## as in Excel ##
  if(model=="highriskAtStart"){
    testcostgain <- (1-prop_highrisk)*numRuledOut.new[[4]]*pwaycost.old[4]
                      + ifelse(cat3TB, 0, (1-prop_highrisk)*numRuledOut.new[[3]]*pwaycost.old[3])
    timegain <- (1-prop_highrisk)*(numRuledOut.new[[4]]*daysToDiag.old[4] - NumDosanjh[4]*ruleouttime)
                      + ifelse(cat3TB, 0, (1-prop_highrisk)*(numRuledOut.new[[3]]*daysToDaig.old[3] - NumDosanjh[3]*ruleouttime))
    testcostloss <- (1-prop_highrisk)*testcost*npatients
    timeloss <- (1-prop_highrisk)*((numRuledOut.new[[1]]+numRuledOut.new[[2]])*FNtime + sum(NumDosanjh[1:2])*ruleouttime) +
                  ifelse(cat3TB, (1-prop_highrisk)*(numRuledOut.new[[3]]*FNtime + NumDosanjh[3]*ruleouttime), 0)
    calcCruleout_hat <- function(totalcostgain, timecostloss, npatients) (totalcostgain-timecostloss)/((1-prop_highrisk)*npatients)
    }else{
    testcostgain <- (1-prop_highrisk)*(1-cat4propfollowup)*
                                  (numRuledOut.new[[4]]*pwaycost.old[4] + ifelse(cat3TB,0, numRuledOut.new[[3]]*pwaycost.old[3]))
    timegain <- numRuledOut.new[[4]]*(1-prop_highrisk)*((1-cat4propfollowup)*daysToDiag.old[4] - (cat4propfollowup*FNtime)) -
                (NumDosanjh[4]*ruleouttime) +
                ifelse(cat3TB,0,
                       numRuledOut.new[[3]]*(1-prop_highrisk)*((1-cat4propfollowup)*daysToDiag.old[3] - (cat4propfollowup*FNtime)) - (NumDosanjh[3]*ruleouttime))
    testcostloss <- testcost*npatients
    timeloss <- (numRuledOut.new[[1]]+numRuledOut.new[[2]])*FNtime*(1-prop_highrisk) + (sum(NumDosanjh[1:2])*ruleouttime) +
                ifelse(cat3TB, numRuledOut.new[[3]]*FNtime*(1-prop_highrisk) + (NumDosanjh[3]*ruleouttime), 0)
    calcCruleout_hat <- function(totalcostgain, timecostloss, npatients) (totalcostgain-timecostloss)/npatients
    }

  timecostgain  <- A*qaly*timegain
  timecostloss  <- A*qaly*timeloss
  totalcostgain <- testcostgain + timecostgain
  totalcostloss <- testcostloss + timecostloss
  totalcostdiff <- totalcostgain - totalcostloss
  Cruleout_hat <- calcCruleout_hat(totalcostgain, timecostloss, npatients)
  INMB <- totalcostdiff/npatients
  INHB <- qaly*(timegain-timeloss)/(365*npatients) - (testcostloss-testcostgain)/(365*A*npatients)
  ICER <- (testcostgain-testcostloss)/(qaly*(timegain-timeloss))


  out.combined <- data.frame(testcost=testcost,
                             sensitivity=sensitivity,
                             specificity=specificity,
                             testcostgain,
                             timegain,
                             timecostgain,
                             testcostloss,
                             timeloss,
                             timecostloss,
                             totalcostgain,
                             totalcostloss,
                             totalcostdiff,
                             INMB,
                             INHB,
                             Cruleout_hat,
                             ICER)


  ##TODO##
  out.costbyDosnajh <- round(data.frame(#scenarioNum=1:length(sensitivity),
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
                                          #pwaycost.old1=pwaycost.old[1],
                                          #pwaycost.old2=pwaycost.old[2],
                                          #pwaycost.old3=pwaycost.old[3],
                                          #pwaycost.old4=pwaycost.old[4],
                                          #cost.new1=cost.new.notsampled[[1]],
                                          #cost.new2=cost.new.notsampled[[2]],
                                          #cost.new3=cost.new.notsampled[[3]],
                                          cost.new4=cost.new.notsampled[[4]],
                                          cost.new.total=cost.new.notsampled.total,
                                          #diff_cost1=diff_cost.old.new.notsampled[[1]],
                                          #diff_cost2=diff_cost.old.new.notsampled[[2]],
                                          #diff_cost3=diff_cost.old.new.notsampled[[3]],
                                          diff_cost4=diff_cost.old.new.notsampled[[4]],
                                          diff_cost.total=diff_cost.old.new.notsampled.total), 2)

  out.diagtimebyDosnajh <- round(data.frame(#scenarioNum=1:length(sensitivity),
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

  # write.csv(out.costbyDosnajh, "../../../output_data/ruleout-costtable-notsampled.csv")
  # write.csv(out.diagtimebyDosnajh, "../../../output_data/ruleout-diagtimetable-notsampled.csv")

  list(out.combined, out.costbyDosnajh, out.diagtimebyDosnajh)
}





## bootstrap sampled patients
sample.IDs <- function(cat){
  sapply(numNotRuledOut.new[[cat]], function(x)
    sample(data$PatientStudyID[data$DosanjhGrouped==cat], size=x, replace = FALSE), simplify = F)
}






calc.ruleOutDiag <- function(cat, N, FNtime, ruleouttime){
  function(nruledin){
    ifelse(cat==4,
           return(ruleouttime*N),
           return((ruleouttime*N) + ((N-nruledin)*FNtime))
    )}
}

calc.daysToDiag <- function(cat){

  N <- NumDosanjh[cat]
  ruleOutDiag <- calc.ruleOutDiag(cat, N, FNtime=0, ruleouttime=1)

  sapply(sampledIDs[[cat]], function(x){
    ((mean(data$start.to.diag[data$PatientStudyID%in%x], na.rm=T)*length(x)) + ruleOutDiag(length(x)))/N
  })
}

calc.ruleOutCost <- function(cat, N, FNcost){
  function(ruleoutcost, nruledin){
    ifelse(cat==4,
           return(ruleoutcost*N),
           return((ruleoutcost*N) + ((N-nruledin)*FNcost))
    )}
}

calc.meancost.new <- function(cat){

  ruleOutCost <- calc.ruleOutCost(cat, N=NumDosanjh[cat], FNcost=0)
  mapply(function(x,y) ((mean(data$totalcost[data$PatientStudyID%in%x], na.rm=T)*length(x)) + ruleOutCost(y, length(x)))/NumDosanjh[cat],
         sampledIDs[[cat]], testcost)
}


##############
## plotting ##
##############

plot.surface_smooth <- function(out){

  axis <- 1:(length(out[[1]]$Cruleout_hat)/length(Ctest))  # unique sens/spec
  d <- data.frame(Sensitivity=out[[1]]$sensitivity[axis], Specificity=out[[1]]$specificity[axis], C=out[[1]]$Cruleout_hat[axis])
  x <- matrix(out[[1]]$Cruleout_hat[axis], ncol=length(unique(out[[1]]$sensitivity[axis])), nrow=length(unique(out[[1]]$sensitivity[axis])))

  ## lattice pkg
  #   image(x=rev(thresh), y=rev(thresh), z=x)
  #   contour(x=rev(thresh), y=rev(thresh), z=x[], add = T)
  # lattice::levelplot(x)
  # lattice::contourplot(x,)

  # ggplot(d, aes(sens, spec, fill=C)) + geom_tile()+scale_fill_gradient(limits = c(100, 500), low = "yellow", high = "red")
  ggplot(d, aes(Sensitivity,Specificity,z=C)) + geom_tile(aes(fill=C))+
    # stat_contour(bins=6, aes(Sensitivity,Specificity,z=C), color="black", size=0.6)+
    stat_contour(aes(Sensitivity,Specificity,z=C, colour="..level.."), color="black", size=0.6, breaks=seq(-100,300,by=10))+
    scale_fill_gradientn(colours=brewer.pal(6,"YlOrRd")) + theme_bw() #, limits=c(40,100)
  # direct.label(v)
}

plot.surface_solid <- function(out){
  # http://stackoverflow.com/questions/10981324/ggplot2-heatmap-with-colors-for-ranged-values

  axis <- 1:(length(out[[1]]$Cruleout_hat)/length(Ctest))  # unique sens/spec
  mat <- data.frame(Sensitivity=out[[1]]$sensitivity[axis], Specificity=out[[1]]$specificity[axis], C=out[[1]]$Cruleout_hat[axis])
  len <- 8
  brks <- cut(mat$C,breaks=seq(0,200,len=len))
  brks <- gsub(","," - ",brks,fixed=TRUE)
  mat$brks <- gsub("\\(|\\]","",brks)  # reformat guide labels

  ggplot(mat,aes(Sensitivity, Specificity))+
    geom_tile(aes(fill=brks))+
    scale_fill_manual("Z", values=brewer.pal(len,"YlOrRd"))+
    # scale_fill_manual(values = colorRampPalette(c("orange", "yellow"))(len))+   #if len>9
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_fixed()
}
