#' Plot Cost-effectiveness Acceptability Curves
#'
#' x-axis of rule-out test sensitivity and grouped
#' by specificity and cost of test
#'
#' @param dat.temp
#'
#' @return NULL
#' @export
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



#' Print a Table of Cost-Effectiveness Acceptability Curve Estimates
#'
#' @param dat.subset
#' @param cost
#' @param sens
#'
#' @return none
#' @export
maketable.CEAC_estimates <- function(dat.subset, cost, sens){

  x <- dat.subset[dat.subset$testcost==cost & dat.subset$sensitivity==sens,]
  print(paste("CEAC>0.5, cost=", cost, ", sens=", sens, " spec=", max(x[x$wmn<=0.5, "specificity"]), sep=""))
  print(paste("CEAC>0.9, cost=", cost, ", sens=", sens, " spec=", max(x[x$wmn<=0.9, "specificity"]), sep=""))
  print(paste("CEAC>0.99, cost=", cost, ", sens=", sens, " spec=", max(x[x$wmn<=0.99, "specificity"]), sep=""))
}


unif_CEAC <- function(x) sum(x>0)/length(x)
beta_CEAC <- function(x) {} ##TODO##



#' Plot Cost-Effectiveness Plane With Contours
#'
#' @param data IDEA study dataset
#' @param N Number of random draws sample points
#' @param Ctest Unit cost of rule-out test
#' @param thresh Sensitivity and specificity of rule-out test (numeric)
#' @param table
#'
#' @return none
#' @examples
#'
#' data("TBdata_clinical_cleaned")
#' plotCEplane(data, N=1000, 400, 0.8)

plotCEplane <- function(data, N = 1000, Ctest, thresh, table=FALSE){

  require(ggplot2)

  npatients <- nrow(data)
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


  if(table==FALSE){
    yearindays <- 365
    return(
      ggplot(data=df, aes(x,y)) + geom_point(colour="grey") + geom_density2d(aes(colour=..level..)) +
        scale_colour_gradient(high="black", low="black") + theme_bw() + theme(legend.position="none") +
        xlab("Health differential") + ylab("Cost differential") +
        ggtitle(paste("Rule-out test cost £", Ctest, "\n and sensitivity=specificity=", thresh, sep="")) +
        ylim(-150,150) + xlim(-5,10) +
        geom_abline(intercept=0, slope=20000/yearindays) +
        geom_abline(intercept=0, slope=30000/yearindays, linetype="dashed")
      )}
  if(table==TRUE){
    return(round(c(
      topright=sum(df$x>0 & df$y>0)/npatients,
      bottomright=sum(df$x>0 & df$y<0)/npatients,
      topleft=sum(df$x<0 & df$y>0)/npatients,
      bottomleft=sum(df$x<0 & df$y<0)/npatients),2))
  }
}
