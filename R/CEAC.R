#' Plot Cost-effectiveness Acceptability Curves
#'
#' x-axis of rule-out test sensitivity and grouped
#' by specificity and cost of test
#'
#' @param dat.temp
#'
#' @return NULL
#'

plotCEAC.xSens.groupSpecCost <- function(dat.temp){
  ## hardcoded y=wmn. could change to mn

  require(ggplot2)
  require(grid)

  ggplot(data=dat.temp,
         aes(x=specificity, y=wmn, group=interaction(sensitivity, testcost), colour=as.factor(sensitivity),
             linetype=as.factor(testcost))) + geom_line(size=2) + theme_bw() +
    ylim(0,1) + xlab("Specificity") + ylab("Probability cost-effective") +
    scale_linetype_discrete(name  ="Rule-out test cost (£)",
                            labels=sort(unique(dat.temp$testcost))) +
    scale_colour_discrete(name  ="Sensitivity",
                          labels=sort(unique(dat.temp$sensitivity))) +
    theme(legend.key.width=unit(3,"line")) +
    geom_hline(aes(yintercept=c(0.9,0.99)), colour="lightgrey")
  #   scale_color_manual("Legend Title\n",labels = c("400", "500"), values = c("black", "grey"))+
  #   theme(legend.title=element_blank())
}



#' Make ruleout table for multiple follow-up times
#'
#' Iterate over a range of follow-up times and row-bind results (potentially a high number of rows).
#'
#' @param FNtime.seq
#' @param Ctest
#'
#' @return out

make.ruleoutTable.FNtime <- function(FNtime.seq=c(7,seq(0,60,by=1)), Ctest){
  out <- NULL
  for (i in FNtime.seq){
    out <- rbind(out,
                 cbind(make.ruleoutTable.pre(FNtime=i, Ctest=Ctest, thresh = seq(from=1, to=0.8, by=-0.005))[[1]], FNtime=i))
  }
  out}


#' Make ruleout table of multiple clinically judged highrisk cut-offs
#'
#' Iterate over a range of cut-off values [0,1] and row-bind results (potentially a high number of rows).
#'
#' @param highrisk.seq
#' @param Ctest
#'
#' @return out

make.ruleoutTable.highrisk <- function(highrisk.seq=seq(0.1,0.99,by=0.01), Ctest){
  out <- NULL
  for (i in highrisk.seq){
    out <- rbind(out,
                 cbind(make.ruleoutTable.pre(prop_highrisk=i, Ctest=Ctest, thresh = seq(from=1, to=0.8, by=-0.005))[[1]], prop_highrisk=i))
  }
  out}


#' Make ruleout table for multiple cohort prevalences
#'
#' Iterate over a range of prevalences 0-100% and row-bind results (potentially a high number of rows).
#'
#' @param cat4percent.seq
#' @param Ctest
#'
#' @return out

make.ruleoutTable.cat4percent <- function(cat4percent.seq=seq(1,100,by=1), Ctest){
  out <- NULL
  for (i in cat4percent.seq){
    out <- rbind(out,
                 cbind(make.ruleoutTable.pre(cat4percent=i, Ctest=Ctest, thresh = seq(from=1, to=0.8, by=-0.005))[[1]], cat4percent=i))
  }
  out}


#' Make ruleout table of bootstrapped samples of patients
#'
#' Iterate over bootstrap samples and row-bind results (potentially a high number of rows).
#'
#' @param nboot
#' @param Ctest
#'
#' @return out

make.ruleoutTable.bootsample <- function(nboot = 2, Ctest){
  out <- NULL
  for (i in 1:nboot){

    bootrows <- sample(1:nrow(data), nrow(data), replace=TRUE)
    bootdata <- data[bootrows,]
    out <- rbind(out,
                 cbind(make.ruleoutTable.pre(data=bootdata, #need to change original function!
                                             Ctest=Ctest, thresh=seq(from=1, to=0.5, by=-0.05))$combinedDosanjh, bootsample=i))
  }
  out}


#' Print a table of Cost-effectiveness Acceptability Curve estimates
#'
#' @param dat.subset
#' @param cost
#' @param sens
#'
#' @return none

maketable.CEAC_estimates <- function(dat.subset, cost, sens){

  x <- dat.subset[dat.subset$testcost==cost & dat.subset$sensitivity==sens,]
  print(paste("CEAC>0.5, cost=", cost, ", sens=", sens, " spec=", max(x[x$wmn<=0.5, "specificity"]), sep=""))
  print(paste("CEAC>0.9, cost=", cost, ", sens=", sens, " spec=", max(x[x$wmn<=0.9, "specificity"]), sep=""))
  print(paste("CEAC>0.99, cost=", cost, ", sens=", sens, " spec=", max(x[x$wmn<=0.99, "specificity"]), sep=""))
}


unif_CEAC <- function(x) sum(x>0)/length(x)
beta_CEAC <- function(x) {} ##TODO##



#' Plot cost-effectiveness plane with contours
#'
#' @param data
#' @param N Number of random draws sample points
#' @param Ctest unit cost of rule-out test
#' @param thresh sensitivity and specificity of rule-out test (numeric)
#'
#' @return none

plotCEplane <- function(data, N, Ctest, thresh){

  require(ggplot2)

  npatients <- nrows(data)
  # grid <- expand.grid(FNtime=rnorm(n = 20, mean=42, sd=10), cat4percent=rbeta(n=3, 5,2), prop_highrisk=rbeta(n=3, 4,6))
  grid <- data.frame(cbind(FNtime=rnorm(n = N, mean=42, sd=10),
                           cat4percent=rbeta(n=N, 5,2)*100,
                           prop_highrisk=rbeta(n=N, 4,6)))
  NCOLgrid <- ncol(grid)
  testrun  <- make.ruleoutTable.pre(data=data, Ctest=50, thresh = 0.99)$combinedDosanjh
  NCOL <- ncol(testrun)
  out  <- matrix(NA, ncol=NCOL, nrow=nrow(grid))
  colnames(out) <- names(testrun)
  out <- data.frame(grid, out)

  for (i in 1:nrow(grid)){

    bootrows <- sample(1:nrow(data), nrow(data), replace=TRUE)
    bootdata <- data[bootrows,]

    out[i,(NCOLgrid+1):NCOL] <- make.ruleoutTable.pre(data=bootdata,
                                                      FNtime=grid$FNtime[i],
                                                      cat4percent=grid$cat4percent[i],
                                                      prop_highrisk=grid$prop_highrisk[i],
                                                      Ctest=Ctest, thresh = thresh)$combinedDosanjh
  }

  df <- data.frame(x = 0.67*(out$timeavoid-out$timeincur)/npatients,
                   y = (out$testcostincur-out$testcostavoid)/npatients)

  ggplot(data=df, aes(x,y)) + geom_point(colour="grey") + geom_density2d(aes(colour=..level..)) +
    scale_colour_gradient(high="black", low="black") + theme_bw() + theme(legend.position="none") +
    xlab("Health differential") + ylab("Cost differential") +
    ggtitle(paste("Rule-out test cost £", Ctest, "\n and sensitivity=specificity=", thresh, sep="")) +
    ylim(-150,150) + xlim(-5,10) +
    geom_abline(intercept=0, slope=55) +
    geom_abline(intercept=0, slope=82, linetype="dashed")
}
