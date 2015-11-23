
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
  points(totalenhanced[baselineRow], 1, pch=8, cex=1.5, lwd=2)

  abline(v=0, lty=3)

  #  ------------------------------------------------------------------------


  # rows <- sample(1:nrow(out$riskprofile$costs), 1000)
  # for (i in rows){
  #   ordered <- order(as.numeric(out$riskprofile$costs[i,]))
  # lines(as.numeric(out$riskprofile$costs[i,])[ordered], cumsum(as.numeric(out$riskprofile$probs[i,])[ordered]), col=rgb(0,0,0,0.1), type="s", pch=16, lwd=5)
  # }

  # lines(max1_
  # lines(max2)
  # lines(max3)
  # lines(max4, type = "p", pch=3)
  #
  # lines(min1)
  # lines(min2)
  # lines(min3)
  # lines(min4, type = "p", pch=3)

  # segments(mean(min1[,1]),min1[1,2],mean(max1[,1]),max1[1,2], lty=2)
  # segments(mean(min2[,1]),min2[1,2],mean(max2[,1]),max2[1,2], lty=2)
  # segments(mean(min3[,1]),min3[1,2],mean(max3[,1]),max3[1,2], lty=2)
  # segments(mean(min4[,1]),min4[1,2],mean(max4[,1]),max4[1,2], lty=2)

  ##

#   mean1 <- aggregate(probs[,1], by=list(costs[,1]), mean)
#   mean2 <- aggregate(probs[,2], by=list(costs[,2]), mean)
#   mean3 <- aggregate(probs[,3], by=list(costs[,3]), mean)
#   mean4 <- aggregate(probs[,4], by=list(costs[,4]), mean)
#
#   sd1 <- aggregate(probs[,1], by=list(costs[,1]), sd)
#   sd2 <- aggregate(probs[,2], by=list(costs[,2]), sd)
#   sd3 <- aggregate(probs[,3], by=list(costs[,3]), sd)
#   sd4 <- aggregate(probs[,4], by=list(costs[,4]), sd)
#
#   upper1 <- lower1 <- mean1
#   upper2 <- lower2 <- mean2
#   upper3 <- lower3 <- mean3
#   upper4 <- lower4 <- mean4
#
#   upper1[,2] <- upper1[,2] + sd1[,2]
#   upper2[,2] <- upper2[,2] + sd2[,2]
#   upper3[,2] <- upper3[,2] + sd3[,2]
#   upper4[,2] <- upper4[,2] + sd4[,2]
#
#   lines(mean1, col="red", lwd=3)
#   lines(mean2, col="red", lwd=3)
#   lines(mean3, col="red", lwd=3)
#   lines(mean4, col="red", lwd=3, type = "p", pch=3)
#
#   lines(upper1, col="red")
#   lines(upper2, col="red")
#   lines(upper3, col="red")
#   lines(upper4, col="red", type = "p", pch=3)
#
#   lower1[,2] <- lower1[,2] - sd1[,2]
#   lower2[,2] <- lower2[,2] - sd2[,2]
#   lower3[,2] <- lower3[,2] - sd3[,2]
#   lower4[,2] <- lower4[,2] - sd4[,2]
#
#   lines(lower1, col="red")
#   lines(lower2, col="red")
#   lines(lower3, col="red")
#   lines(lower4, col="red", type = "p", pch=3)
#
#   segments(mean(lower1[,1]),lower1[1,2],mean(upper1[,1]),upper1[1,2], col="red")
#   segments(mean(lower2[,1]),lower2[1,2],mean(upper2[,1]),upper2[1,2], col="red")
#   segments(mean(lower3[,1]),lower3[1,2],mean(upper3[,1]),upper3[1,2], col="red")
#   segments(mean(lower4[,1]),lower4[1,2],mean(upper4[,1]),upper4[1,2], col="red")
#
#   segments(mean(lower3[,1]),mean3[1,2],mean(lower4[,1]),mean4[1,2])
#   segments(mean(lower2[,1]),mean2[1,2],mean(lower4[,1]),mean4[1,2])
#   segments(mean(lower1[,1]),mean1[1,2],mean(lower2[,1]),mean2[1,2])





}

