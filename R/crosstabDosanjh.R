#' Create contingency table by Dosanjh categories of data and ideal pathways
#'
#' \code{crosstabDosanjh} returns
#'
#' @param out output from \code{createDecisionTree}
#' @return table

crosstabDosanjh <- function(out){

  # addmargins(table(out$t1$data$treat1, out$t1$data$treat1.est))
  # chisq.test(ftable(out$t1$data$TBconfirmed, out$t1$data$treat1))
  # ftable(out$t1$data$TBconfirmed, out$t1$data$treat1, out$t1$data$treat1.est)
  xtab <- addmargins(ftable(out$t1$data$TBconfirmed, out$t1$data$treat1, out$t1$data$treat1.est))
  write.csv(xtab, file="./treat1crosstab.csv")

  # ftable(out$treat1obs_Sxobs$data$TBconfirmed, out$treat1obs_Sxobs$data$treat2.ideal, out$treat1obs_Sxobs$data$treatResponse)
  xtab.treat2 <- addmargins(ftable(out$treat1est_Sxest$data$TBconfirmed, out$treat1est_Sxest$data$treatResponse,
                                   out$treat1est_Sxest$data$treat2.ideal))
  write.csv(xtab.treat2, file="./treat2crosstab.csv")

  # table(out$treat1est_Sxest$data$DosanjhGrouped, out$treat1obs_Sxest$data$DosanjhIdeal)
  xtab.Dosanjh.est <- addmargins(ftable(out$treat1est_Sxest$data$TBconfirmed, out$treat1est_Sxest$data$DosanjhGrouped,
                                    out$treat1est_Sxest$data$DosanjhIdeal))
  write.csv(xtab.Dosanjh.est, file="./Dosanjhcrosstab_est.csv")

  xtab.Dosanjh <- addmargins(ftable(out$treat1obs_Sxest$data$TBconfirmed, out$treat1obs_Sxest$data$DosanjhGrouped,
                                    out$treat1obs_Sxest$data$DosanjhIdeal))
  write.csv(xtab.Dosanjh, file="./Dosanjhcrosstab.csv")
}


