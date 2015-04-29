#' plot.FollowUpChart
#'
#' \code{plot.FollowUpChart} returns joined tables of inputs and predicted outcomes
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

plot.FollowUpChart <- function(data, col=1, lty=1, lwd=2, xlim, xlab="Time (days)", ylab="", order="bydrugstart",
                               legend=TRUE, dropNA=TRUE, axisID=TRUE, ...){
    ## plot representation of patient journeys
    ##
    ## dropNA: should records with NAs for test to drug times be included?
    ## axisID: include patient ID code as y-axis label or not

    c25 <- c("dodgerblue2","#E31A1C", # red
             "green4",
             "#6A3D9A", # purple
             "#FF7F00", # orange
             "black","gold1",
             "skyblue2","#FB9A99", # lt pink
             "palegreen2",
             "#CAB2D6", # lt purple
             "#FDBF6F", # lt orange
             "gray70", "khaki2",
             "maroon","orchid1","deeppink1","blue1","steelblue4",
             "darkturquoise","green1","yellow4","yellow3",
             "darkorange4","brown")

  if (missing(xlim)) xlim <- c(min(data),max(data))

    if(order=="bytreatment"){
        data <- data[order(data$TBDrug_diff, data$testDrug_diff, data$testDiagCon_diff),]
    }else if(order=="bydrugstart"){
        data <- data[order(data$testDrug_diff, data$testDiagCon_diff),]}

    if(dropNA){N <- nrow(data[!is.na(data$testDrug_diff),])
    }else{N <- nrow(data)}

  ## empty plot
  if(axisID==FALSE){
      plot(0,0,type="n",ylim=c(1,N+1),xlim=xlim,xlab=xlab,ylab="", bg="white", ...)
  }else{plot(0,0,type="n",ylim=c(1,N+1),xlim=xlim,xlab=xlab,ylab="", bg="white", yaxt="n", las=1, cex.axis=0.5, ...)
        axis(2, at=1:N, labels=data$PatientStudyID[1:N], las=1, cex.axis=0.5)}

  nix <- lapply(1:N, function(i){
                        segments(0, i, data$testDrug_diff[i], i, col='lightgray',lty=lty,lwd=lwd)
                        segments(data$testDrug_diff[i],i, data$testDrug_diff[i]+data$TBDrug_diff[i],i ,col=col,lty=lty,lwd=lwd)
  })

  circle <- rep.int(19, times = N)
  circle[data$TBcultCens] <- 21

    cultCols <- rep("red",N)
    cultCols[data$TBcult[1:N]=="POSITIVE"] <- "darkgreen"
    cultCols[is.na(data$TBcult[1:N])] <- "blue"
    cultCols[data$TBcult[1:N]=="Not taken"] <- "orange"
    cultCols[data$TBcult[1:N]=="INDETERMINATE"] <- "yellow"

    ticktime <- apply(cbind(data$testDrug_diff[1:N]+60, rep(60,N)), 1, max, na.rm=TRUE)

    points(ticktime, 1:N, pch=15, cex=1, col=c25[as.numeric(as.factor(data$DiagoutcomeGrouped))])    #[1:N]+2)  #tick-marks
    points(data$testDrug_diff[1:N]+183, 1:N, pch=4, cex=1, col="black")
    points(data$testDrug_diff[1:N]+365, 1:N, pch=4, cex=1, col="black")
    points(rep(0,N), 1:N, pch='l', col="lightgrey")
    points(data$testCult_diff[1:N], 1:N, pch=circle, cex=1, col=cultCols)# col="red")

  if (legend==TRUE){
    Diagoutcome.levels <- levels(as.factor(data$DiagoutcomeGrouped))
    legend("bottomright",
           c(Diagoutcome.levels,"First test to treatment start","Treatment period","Culture report: negative","Culture report: positive","6 and 12 months treatment duration"),
           col=c(c25[1:length(Diagoutcome.levels)],
                #c(1:length(Diagoutcome.levels)+2,
                 'lightgray',1,'red','darkgreen',1),
           pch=c(rep(15,length(Diagoutcome.levels)),3,3,19,19,4))}
}

