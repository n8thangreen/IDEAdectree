
#' Incremental Net Monetary Benefit for IDEA enhanced diagnostic pathway
#'
#' This a kind of super slimmed-down version of \code{make.ruleoutTable.pre()}.
#'
#' @param spec
#' @param sens
#' @param pathcost Routine in-practice total cost (time and tests) for a non-TB patient (433+11*55*0.67)
#' @param Ctest Unit cost of rule-out test
#' @param A Cost-effectiveness threshold (per year)
#' @param q QALY
#' @param prev Cohort prevalence. An 'uncertain' parameter
#' @param FNtime An 'uncertain' parameter (days)
#' @param prop_highrisk Threshold high risk predictive probability. An 'uncertain' parameter
#'
#' @return INMB

INMB <- function(spec=0.9, sens=0.90, Ctest=100, pathcost=838, A=30000, q=0.67,
                 prev, FNtime, prop_highrisk){

  # require(data.table)

  stopifnot(prev>=0, prev<=1)
  stopifnot(prop_highrisk>=0, prop_highrisk<=1)
  stopifnot(FNtime>=0)

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


#' Calculate Expected INMB
#'
#' @param spec
#' @param sens
#' @param Ctest
#'
#' @return INMB sample

calcExpectedINMB <- function(spec=0.9, sens=0.9, Ctest=200){

  grid <- data.frame(cbind(FNtime=rnorm(n=N, mean=42, sd=10),
                           cat4percent=rbeta(n=N, 5,2)*100,
                           prop_highrisk=rbeta(n=N, 4,6)))

  out <- NA
  for (i in 1:nrow(grid)){
    out[i] <- INMB(spec, sens, Ctest,
                   FNtime=grid$FNtime[i], prev=1-(grid$cat4percent[i]/100), prop_highrisk=grid$prop_highrisk[i])
  }
  out
}


#' Make Expected INMB Table
#'
#' @return matrix
# write.csv(makeExpectedINMBtable(), "temp.csv")
makeExpectedINMBtable <- function(){

  input <- data.frame(cbind(spec=rep(c(0.8,0.9,0.99), each=3),
                            sens=rep(c(0.8,0.9,0.99), each=3),
                            Ctest=c(400,300,200)))
  out <- outmed <- out1 <- out3 <- NULL
  for (i in 1:nrow(input)){

    summaryINMB <- summary(calcExpectedINMB(spec=input$spec[i], sens=input$sens[i], Ctest=input$Ctest[i]))
    # out <- c(out, paste(summaryINMB["Median"], " [", summaryINMB["1st Qu."], ",", summaryINMB["3rd Qu."], "]", sep=""))
    outmed <- c(outmed, summaryINMB["Median"])
    out1 <- c(out1, summaryINMB["1st Qu."])
    out3 <- c(out3, summaryINMB["3rd Qu."])
  }
  # return(rbind(t(input), out))
  return(rbind(t(input), outmed, out1, out3))
}


#' Calculate Expected Value of Information
#'
#' @return array
#'
#' @seealso \code{\link{plot.EVPI}}

calc.EVPI <- function(N=5000){

  ## expected value over all unknown parameters

  seqSPEC  <- seq(0.8, 1, by=0.02)
  seqSENS  <- c(0.9, 0.99)
  seqCtest <- c(550, 600)

  lseqSPEC <- length(seqSPEC)
  EVunknown <- rep(NA, lseqSPEC)
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
#' @seealso \code{\link{calc.EVPI}}


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
}


