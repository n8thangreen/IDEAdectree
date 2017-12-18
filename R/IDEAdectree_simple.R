
#' IDEA Decision Tree Calculation
#'
#' This is a parred down version of the full model that is used in the paper.
#' This is specifically for the HTA chapter.
#' It probabilistically incorporates sampling variability (and standard pathway cost and time to diagnosis).
#'
#' @param data IDEA study data
#' @param nsim Number of sample points. Default: 1000
#' @param costDistns List of distribution names and parameter values for each test/procedure
#' @param prevalence As a probability. Default: 0.25
#' @param cutoff Clinical judgement threshold
#' @param FNtime False negative follow-up time
#' @param FNdist Should false negative time to follow-up distribution be used (logical)
#' @param SENS Sensitivity Rule-out test
#' @param SPEC Specificity Rule-out test
#' @param SENSvar Sensitivity variance
#' @param SPECvar Specificiaty variance
#' @param c.ruleout Rule-out test unit cost
#' @param name.ruleout Name of rule-out test to get at distribution
#' @param quant Quantile value of time to diagnosis and costs
#' @param WTP Willingness-to-pay
#' @param utility QALY adjustment utility due to active TB
#' @param N Number of patients. The number in the data is used as default.
#' @param wholecohortstats Should the output stats be the total or per patient
#'
#' @return Health and cost realisations
#'
#' @examples
#'
#' library(bcea)
#' dat1 <- IDEAdectree.simple(data, cutoff = 0.4, specificity = 0.8)
#' dat2 <- IDEAdectree.simple(data, cutoff = 0.4, specificity = 0.9)
#' dat3 <- IDEAdectree.simple(data, cutoff = 0.4, specificity = 0.99)
#'
#' dat$e <- cbind(dat1$e, dat2$e[,2], dat3$e[,2])
#' dat$c <- cbind(dat1$c, dat2$c[,2], dat3$c[,2])
#'
#' intlabels <- c("Current", "Enhanced specificity=0.8", "Enhanced specificity=0.9", "Enhanced specificity=0.99")
#'
#' m <- bcea(e=dat$e, c=dat$c, ref=1, interventions = intlabels)
#' contour2(m, wtp=WTP, graph = "ggplot2", ICER.size=2, pos=c(0.9,0.1), xlim=c(-5,20), ylim=c(-100,500)) + ggtitle("")
#' summary(m)

IDEAdectree.simple <- function(data,
                               nsim = 1000,
                               costDistns = COST.distns.allerror,  #COST.distns
                               prevalence = 0.25,
                               cutoff = 1,
                               FNtime = 42,
                               FNdist = TRUE,
                               SENS = 0.9,
                               SPEC = 0.9,
                               SENSvar = 0.005,
                               SPECvar = 0.005,
                               c.ruleout = 100,
                               name.ruleout = NA,
                               quant = 0.5,
                               WTP = 20000, #30000
                               utility = NA,
                               N = nrow(data),
                               wholecohortstats = FALSE) {
  # require(assertive)
  require(triangle)

  stopifnot(name.ruleout %in% c(NA, names(costDistns)))
  stopifnot(nsim > 0,
            prevalence >= 0,
            prevalence <= 1,
            cutoff >= 0,
            cutoff <= 1,
            FNtime >= 0,
            SENS >= 0,
            SENS <= 1,
            SPEC >= 0,
            SPEC <= 1,
            c.ruleout >= 0,
            WTP >= 0)
  stopifnot(quant >= 0, quant <= 1)

  A <- WTP/365.25  #per day
  #QALY adjustment absolute utility
  e <- c <- NULL

  FNdens <-
    as.numeric(data$start.to.FU) %>%
    subset(data$VisitFU == "2 month FU" &
             data$DosanjhGrouped %in% c(1,2) &
             !is.na(data$start.to.FU) &
             data$start.to.FU <= 200 &
             data$start.to.FU > 0) %>%
    density(from = 0, bw = 10)

  Fx <- cumsum(FNdens$y)/sum(FNdens$y)

  sens.betaparams <- MoM.beta(xbar = SENS, vbar = SENSvar)
  spec.betaparams <- MoM.beta(xbar = SPEC, vbar = SPECvar)
  sensitivity <- rbeta(n = nsim, shape1 = sens.betaparams$a, shape2 = sens.betaparams$b)
  specificity <- rbeta(n = nsim, shape1 = spec.betaparams$a, shape2 = spec.betaparams$b)

  visit1cost <- rgamma(n = nsim, shape = 53.3, scale = 4.52)
  visit2cost <- rgamma(n = nsim, shape = 18.78, scale = 7.62)

  if (is.na(utility)) {
    utility <- rtriangle(n = nsim, a = 0.69, b = 0.89)  #1-QALY loss i.e. relative to non-TB (<1)
  }else{
    utility <- rep(utility, time = nsim)}

  ## don't include generic tests costs
  costDistns$PET$params["mean"] <- 0
  costDistns$MRI$params["mean"] <- 0
  costDistns$CT$params["mean"]  <- 0


  for (i in 1:nsim) {

    rcosts <- sample.distributions(costDistns)
    totalcost <- calcPatientCostofTests(data, COSTS = rcosts)

    if ("PatientWgt" %in% names(data)) {
      weight <- mean(data$PatientWgt, na.rm = TRUE)
    }else{
      weight <- 67.98}  #kg

    whoCat4Treated <- !is.na(data$TBDrugStart.min) & data$DosanjhGrouped == "4"

    treatment.days <- 60
    twomonthTreatCost <- treatment.days * ((rifamicin.cost_qty*rifamicin.mg_day)/(rifamicin.pill_qty*rifamicin.mg_pill) +
                                             (isoniazid.cost_qty*isoniazid.mg_day)/(isoniazid.pill_qty*isoniazid.mg_pill) +
                                             (pyrazinamide.cost_qty*pyrazinamide.mg_day)/(pyrazinamide.pill_qty*pyrazinamide.mg_pill) +
                                             (ethambutol.cost_qty*ethambutol.mg_day_kg*weight)/(ethambutol.pill_qty*ethambutol.mg_pill))

    totalcost[whoCat4Treated] <- totalcost[whoCat4Treated] + twomonthTreatCost

    if (!is.na(name.ruleout)) c.ruleout <- rcosts[[name.ruleout]]

    # param.distns <- list(t.ruleout=list(distn="unif", params=c(min=2, max=7)))
    # t.ruleout <- sample.distributions(param.distns)
    t.ruleout <- runif(1, min = 2, max = 14)
    h.ruleout <- utility[i] * t.ruleout

    if (FNdist) {h.FN <- utility[i] * FNdens$x[sum(runif(1) > Fx)]
    }else{h.FN <- utility[i]*FNtime}

    TB  <- rbinom(n = 1, size = N, prob = prevalence)
    nTB <- N - TB

    TBhighrisk <- rbinom(n = 1, size = TB, prob = 1 - pbeta(cutoff, 7, 3))
    TBlowrisk  <- TB - TBhighrisk

    nTBhighrisk <- rbinom(n = 1, size = nTB, prob = 1 - pbeta(cutoff,3,7))
    nTBlowrisk  <- nTB - nTBhighrisk

    TBpos <- rbinom(n = 1, size = TBlowrisk, prob = sensitivity[i])
    TBneg <- TBlowrisk - TBpos

    nTBpos <- rbinom(n = 1, size = nTBlowrisk, prob = 1 - specificity[i])
    nTBneg <- nTBlowrisk - nTBpos

    ## final subpopulation sizes
    pop  <- c(TBhighrisk, nTBhighrisk, TBpos, TBneg, nTBpos, nTBneg)
    stopifnot(sum(pop) == N)

    ## current time and cost estimates ##

    sboot.nonTB <- sample(which(data$DosanjhGrouped == 4), replace = TRUE)
    sboot.TB <- sample(which(data$DosanjhGrouped %in% c(1,2,3)), replace = TRUE)

    c.std.nonTB <- quantile(totalcost[sboot.nonTB], probs = quant, na.rm = TRUE)
    c.std.TB <- quantile(totalcost[sboot.TB], probs = quant, na.rm = TRUE)
    start.to.diag.nonTB <- quantile(data$start.to.diag[sboot.nonTB], probs = quant, na.rm = TRUE)
    start.to.diag.TB <- quantile(data$start.to.diag[sboot.TB], probs = quant, na.rm = TRUE)

    # h.std <- utility[i]*start.to.diag
    h.std.TB <- utility[i]*start.to.diag.TB
    h.std.nonTB <- utility[i]*start.to.diag.nonTB

    ## outcomes
    # cost <- c(c.std, c.std, c.std+c.ruleout, c.std+c.ruleout, c.std+c.ruleout, c.ruleout)
    # health <- c(h.std, h.std, h.std+h.ruleout, h.std+h.ruleout+h.FN, h.std+h.ruleout, h.ruleout)

    cost   <- visit1cost[i] + c(c.std.TB, c.std.nonTB, c.std.TB + c.ruleout, c.std.TB + c.ruleout + visit2cost[i], c.std.nonTB + c.ruleout, c.ruleout)
    health <- c(h.std.TB, h.std.nonTB, h.std.TB + h.ruleout, h.std.TB + h.ruleout + h.FN, h.std.nonTB + h.ruleout, h.ruleout)

    Ec.old <- (visit1cost[i]*N) + (c.std.TB*TB + c.std.nonTB*nTB)     #THIS IS A BIT OF A PROBLEM FOR VARYING PREVALENCE IN ORDER TO BE A COMPARISON FOR ALL...
    Ee.old <- (h.std.TB*TB) + (h.std.nonTB*nTB)                       #SIMILARLY THE SAMPLED CURRENT TIMES AND COSTS

    if (wholecohortstats) N <- 1

    Ec.old <- Ec.old/N
    Ee.old <- Ee.old/N


    ## expected values
    if (length(pop) == length(health)) {  ##TODO## fix this bug. quick fix
      if (!is.na(Ec.old) & !is.na(Ee.old)) {
        e <- rbind(e, c(Ee.old, (pop %*% health)/N))
        c <- rbind(c, c(Ec.old, (pop %*% cost)/N))}
    }else{
      e <- rbind(e, e[nrow(e), ])
      c <- rbind(c, c[nrow(c), ])
    }
  }

  return(list(e = e,
              c = c))
}


