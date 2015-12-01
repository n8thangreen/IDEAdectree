#
# cat4percent <- seq(0, 100, length = 11)
# cat4percentPROB <- dbeta(cat4percent/100, 5,2)
# cat4percentPROB <- cat4percentPROB/sum(cat4percentPROB)
# thresh <- seq(0, 1, 0.01)
# INMBoptimal <- EPVI <- mat <- NULL
# resmax <- 0
#
# for (i in cat4percent){
#
#   out <- make.ruleoutTable.pre(cat4percent = i, thresh = thresh, Ctest=400)$combinedDosanjh
#
#   INMBenhanced <- out$INMB[out$sensitivity==0.8 & out$specificity%in%thresh]
#   mat <- rbind(mat, INMBenhanced)
# }
#
# INMB.EV <- pmax(0, cat4percentPROB%*%mat)
# INMB.EVwithPI <- cat4percentPROB%*%matrix(pmax(0, mat), nrow=nrow(mat))
# EVPI <- INMB.EVwithPI - INMB.EV
#
# plot(thresh, EVPI, xlab="specificity", type="o", xlim=c(0.7,1))

# https://en.wikipedia.org/wiki/Expected_value_of_perfect_information
# which (single) decision would maximise the expected total utility?
# what decision for each cat4percent would maximise the utility?

##TODO##
# do the same for ctest on xaxis too...






# using oakley paper

#  ------------------------------------------------------------------------

library(data.table)


#' Incremental Net Monetary Benefit for IDEA enhanced diagnostic pathway
#'
#' @param spec
#' @param sens
#' @param pathcost Routine in-practice total cost (time and tests) for a non-TB patient (433+11*55*0.67)
#' @param Ctest
#' @param A Cost-effectiveness threshold (per year)
#' @param q
#' @param prev
#' @param FNtime
#' @param prop_highrisk
#'
#' @return INMB

INMB <- function(spec=0.9, sens=0.90, pathcost=838, Ctest=100, A=20000, q=0.67,
                 prev, FNtime, prop_highrisk){

  A <- A/365
  ruleoutcost <- Ctest + A*q

  ECDF.TB <- stats::ecdf(data$riskfacScore[data$DosanjhGrouped%in%c(1,2)])
  ECDF.nonTB <- stats::ecdf(data$riskfacScore[data$DosanjhGrouped==4])

  spec.clinical <- ECDF.nonTB(prop_highrisk)
  sens.clinical <- 1-ECDF.TB(prop_highrisk)

  return((1-prev)*spec*spec.clinical*pathcost -
           prev*(1-sens.clinical)*(ruleoutcost+(1-sens)*A*q*FNtime) -
           (1-prev)*spec.clinical*ruleoutcost)
}


#
# ## expected value over all unknown parameters
# N <- 40
# seqSPEC <- seq(0.8, 1, by=0.02)
# EVunknown <- rep(NA, length(seqSPEC))
# out <- rep(NA, N)
#
# for (j in 1:length(seqSPEC)){
#
#   grid <- data.frame(cbind(FNtime = rnorm(n = N, mean=42, sd=1),
#                            prev = 1-rbeta(n=N, 5,2),
#                            prop_highrisk = rbeta(n=N, 4,6)))
#   for (i in 1:N){
#
#     out[i] <- INMB(spec=seqSPEC[j],
#                    FNtime=grid$FNtime[i], prev=grid$prev[i], prop_highrisk=grid$prop_highrisk[i])
#   }
#   EVunknown[j] <- max(0, mean(out))
# }
#
#
#
# ## partial expected value of perfect information
#
# res <- rep(NA, length(seqSPEC))
# EVres <- matrix(NA, nrow = nrow(grid), ncol = length(seqSPEC))
# out <- rep(NA, N)
#
# for (k in 1:nrow(grid)){
#
#   for (j in 1:length(seqSPEC)){
#
#     grid <- data.frame(cbind(FNtime = rnorm(n = N, mean=42, sd=1),
#                              prev = 1-rbeta(n=N, 5,2),
#                              prop_highrisk = rbeta(n=N, 4,6)))
#     for (i in 1:N){
#
#       out[i] <- INMB(spec=seqSPEC[j],
#                      FNtime=grid$FNtime[i], prev=grid$prev[k], prop_highrisk=grid$prop_highrisk[i])
#     }
#     res[j] <- max(0, mean(out))
#   }
#   EVres[k,] <- res
# }
# EVPI <- apply(EVres, 2, mean)
#
# plot(seqSPEC, EVPI - EVunknown, type="o", xlab="Specificity")
#
# plot(seqSPEC, EVPI, ylim=c(-5,5), type="l", xlab="Specificity")
# lines(seqSPEC, EVunknown, col="red")



## hard-code parameter expected values
#  ------------------------------------------------------------------------


#' Calculate Expected Value of Information
#'
#' @return array
#'
#' @seealso \code{\link{plot.EVPI}}

calc.EVPI <- function(){

  ## expected value over all unknown parameters

  seqSPEC <- seq(0.8, 1, by=0.02)
  lseqSPEC <- length(seqSPEC)
  EVunknown <- rep(NA, lseqSPEC)
  seqSENS <- c(0.9, 0.99)
  seqCtest <- c(550, 600)
  EVunknown.sens <- NULL

  for (Ctest in seqCtest){
    for (sens in seqSENS){
      for (j in 1:lseqSPEC){

        out <- INMB(spec=seqSPEC[j], sens=sens, Ctest=Ctest,
                    FNtime=42, prev=0.25, prop_highrisk=0.4)
        EVunknown[j] <- max(0, out)
      }
      EVunknown.sens <- rbind(EVunknown.sens, data.frame(Ctest, sens, spec=seqSPEC, EVunknown))
    }
  }

  ## partial expected value of perfect information (EVPI)

  N <- 4000
  res <- rep(NA, lseqSPEC)
  EVres <- matrix(NA, nrow = N, ncol = lseqSPEC)
  EVPIsens <- NULL

  for (Ctest in seqCtest){
    for (sens in seqSENS){
      for (k in 1:N){

        for (j in 1:lseqSPEC){

          prev <- 1 - rbeta(n=1, 5,2)
          out <- INMB(spec=seqSPEC[j], sens=sens, Ctest=Ctest,
                      FNtime=42, prev=prev, prop_highrisk=0.4)
          res[j] <- max(0, out)
        }
        EVres[k,] <- res
      }
      EVPI <- apply(EVres, 2, mean)
      EVPIsens <- rbind(EVPIsens, data.frame(Ctest, sens, spec=seqSPEC, EVPI))
    }
  }

  EVPIsens <- merge(EVPIsens, EVunknown.sens, by=c("sens", "spec", "Ctest"))
  EVPIsens$partialEVPI <- EVPIsens$EVPI - EVPIsens$EVunknown
  EVPIsens
}


#' Plot Expected Value of Information curves
#'
#' @param EVPIsens output from \code{calc.EVPI}
#' @return NULL
#'
#' @seealso \code{\link(calc.EVPI}}


plot.EVPI <- function(EVPIsens){

  require(ggplot2)
  require(grid)

  ggplot(data=EVPIsens, aes(x=spec, y=partialEVPI, group=interaction(sens, Ctest), colour=as.factor(sens),
                            linetype=as.factor(Ctest))) + geom_line() + theme_bw() +
    xlab("Specificity") + ylab("EVPI (£)") +
    scale_linetype_discrete(name  ="Rule-out test cost (£)",
                            labels=sort(unique(EVPIsens$Ctest))) +
    scale_colour_discrete(name  ="Sensitivity",
                          labels=sort(unique(EVPIsens$sens))) +
    theme(legend.key.width=unit(3,"line"))


  #
  # plot(seqSPEC, EVPI - EVunknown, type="o", xlab="Specificity", ylab="EVPI")
  #
  # plot(seqSPEC, EVPI, ylim=c(-5,25), type="l", xlab="Specificity")
  # lines(seqSPEC, EVunknown, col="red")
}


