
#' Plot risk profile discrete distribution
#'
#' Cost and probability of final outcomes (along each distinct path) in decision tree.
#'
#' @param data
#' @param prop_highrisk
#' @param FNtime
#'
#' @return none

plot.riskprofile <- function(data, prop_highrisk, FNtime){

  out <- make.ruleoutTable.pre(prop_highrisk = prop_highrisk, FNtime=FNtime)

  probs <- sapply(out$riskprofile$probs/rowSums(out$riskprofile$probs), as.numeric) ##TODO## check that not summing to 1 is just rounding error
  costs <- sapply(out$riskprofile$costs, as.numeric)

  plot(1, 1,
       ylim=c(0,1), xlim=c(-500,2000), type="n",
       xlab="Difference between current and enhanced diagnostic pathways total cost (£)", ylab="Probability")

  max1 <- aggregate(probs[,1], by=list(costs[,1]), max)
  max2 <- aggregate(probs[,2], by=list(costs[,2]), max)
  max3 <- aggregate(probs[,3], by=list(costs[,3]), max)
  max4 <- aggregate(probs[,4], by=list(costs[,4]), max)

  min1 <- aggregate(probs[,1], by=list(costs[,1]), min)
  min2 <- aggregate(probs[,2], by=list(costs[,2]), min)
  min3 <- aggregate(probs[,3], by=list(costs[,3]), min)
  min4 <- aggregate(probs[,4], by=list(costs[,4]), min)

  polygon(c(max1[,1], rev(min1[,1])), c(max1[,2], min1[,2]), col="grey", border=NA)
  polygon(c(max2[,1], rev(min2[,1])), c(max2[,2], min2[,2]), col="grey", border=NA)
  polygon(c(max3[,1], rev(min3[,1])), c(max3[,2], min3[,2]), col="grey", border=NA)
  polygon(c(max4[,1], rev(min4[,1])), c(max4[,2], min4[,2]), col="grey", border=NA)

  baselineRow <- with(out$combinedDosanjh, sensitivity==0.9 & specificity==0.9 & testcost==200)
  points(costs[baselineRow,], probs[baselineRow,], pch=19)

  segments(costs[baselineRow,3], probs[baselineRow,3], costs[baselineRow,4], probs[baselineRow,4])
  segments(costs[baselineRow,2], probs[baselineRow,2], costs[baselineRow,4], probs[baselineRow,4])
  segments(costs[baselineRow,1], probs[baselineRow,1], costs[baselineRow,2], probs[baselineRow,2])

  totalenhanced <- diag(probs%*%t(costs))
  lines(c(min(totalenhanced), max(totalenhanced)), c(1,1), col="grey")

  points(0,1, pch=8, cex=1.5, lwd=2, col="red")
  points(totalenhanced[baselineRow], 1, pch=0, cex=2, lwd=2)

  abline(v=0, lty=3)
}



