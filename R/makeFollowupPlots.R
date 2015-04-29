#' makeFollowupPlots
#'
#' \code{makeFollowupPlots} returns joined tables of inputs and predicted outcomes
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

makeFollowupPlots <- function(data = "../../output_data/TBdata_prepd_CLEANED_190315.RData"){
  ##
  ##

  load(data)

  ## remove individuals
  datat <- data
  datat <- datat[datat$HIVpos==FALSE,]

  setwd("../../output_plots/data_indiv_pathways")

  xlim.max <- 300
  # x11()

  postscript("Dosanjh4_followup.eps")
  par(mfrow=c(2,2))
  plot.FollowUpChart(datat[datat$Dosanjh=="4A",], xlim=c(-1,xlim.max), main="Dosanjh 4A")
  plot.FollowUpChart(datat[datat$Dosanjh=="4B",], xlim=c(-1,xlim.max), main="Dosanjh 4B")
  plot.FollowUpChart(datat[datat$Dosanjh=="4C",], xlim=c(-1,xlim.max), main="Dosanjh 4C")
  plot.FollowUpChart(datat[datat$Dosanjh=="4D",], xlim=c(-1,xlim.max), main="Dosanjh 4D")
  dev.off()

  postscript("Dosanjh1-3_followup.eps")
  par(mfrow=c(2,2))
  plot.FollowUpChart(datat[datat$Dosanjh=="3",], xlim=c(-1,xlim.max), main="Dosanjh 3", order="bydrugstart")
  plot.FollowUpChart(datat[datat$Dosanjh=="2",], xlim=c(-1,xlim.max), main="Dosanjh 2", order="bydrugstart")
  plot.FollowUpChart(datat[datat$Dosanjh=="1",], xlim=c(-1,xlim.max), main="Dosanjh 1", order="bydrugstart")
  dev.off()

  postscript("PTBorEPTB_followup.eps")
  par(mfrow=c(1,2))
  plot.FollowUpChart(datat[datat$EPTBorPTB=="PTB",], xlim=c(-1,xlim.max), main="PTB")
  plot.FollowUpChart(datat[datat$EPTBorPTB=="EPTB;PTB",], xlim=c(-1,xlim.max), main="EPTB")
  dev.off()

  postscript("HIV_followup.eps")
  par(mfrow=c(1,2))
  plot.FollowUpChart(datat[datat$HIVpos==TRUE,], xlim=c(-1,xlim.max), main="HIV positive")
  plot.FollowUpChart(datat[datat$HIVpos==FALSE,], xlim=c(-1,xlim.max), main="HIV negative")
  dev.off()

  postscript("Homeless_followup.eps")
  par(mfrow=c(1,2))
  plot.FollowUpChart(datat[datat$CurrHomeless==TRUE,], xlim=c(-1,xlim.max), main="Homeless")
  plot.FollowUpChart(datat[datat$CurrHomeless==FALSE,], xlim=c(-1,xlim.max), main="Not homeless")
  dev.off()

  postscript("TBcontact_followup.eps")
  par(mfrow=c(1,2))
  plot.FollowUpChart(datat[datat$TBcont==TRUE,], xlim=c(-1,xlim.max), main="TB contact")
  plot.FollowUpChart(datat[datat$TBcont==FALSE,], xlim=c(-1,xlim.max), main="No TB contact")
  dev.off()


  ethn.lookup <- read.csv("../../raw_data/IDEA_clinical_database/ethnic_lookup.csv")
  for (i in unique(data$Ethn[!is.na(data$Ethn)])){
    names <- paste(ethn.lookup[ethn.lookup$EthnCode==i,"Ethnclass"], ethn.lookup[ethn.lookup$EthnCode==i,"Ethndetail"])
    names <- gsub("/","-", names)
    postscript(paste(names, "followup.eps"))
    plot.FollowUpChart(datat[datat$Ethn==i,], xlim=c(-1,xlim.max), main=names)
    dev.off()
  }


  for (i in levels(datat$WHOcut)){
    postscript(paste(i, "followup.eps"))
    plot.FollowUpChart(datat[datat$WHOcut==i,], xlim=c(-1,xlim.max), main=i)
    dev.off()
  }

  for (i in levels(datat$Region)){
    postscript(paste(i, "followup.eps"))
    plot.FollowUpChart(datat[datat$Region==i,], xlim=c(-1,xlim.max), main=i)
    dev.off()
  }

  for (i in levels(datat$Sitename)){
    postscript(paste(i, "followup.eps"))
    plot.FollowUpChart(datat[datat$Sitename==i,], xlim=c(-1,xlim.max), main=i)
    dev.off()
  }


  ## all Smear positives/negatives
  postscript("Smear_followup.eps")
  par(mfrow=c(1,2))
  plot.FollowUpChart(datat[datat$Smear=="POSITIVE",], xlim=c(-1,300), main="Smear positive", order="bytreatment")
  plot.FollowUpChart(datat[datat$Smear=="NEGATIVE",], xlim=c(-1,400), main="Smear negative")
  dev.off()


  ## all IGRA positives/negatives
  postscript("IGRA_followup.eps")
  par(mfrow=c(1,2))
  plot.FollowUpChart(datat[datat$IGRA=="POSITIVE",], xlim=c(-1,700), main="IGRA positive")
  plot.FollowUpChart(datat[datat$IGRA=="NEGATIVE",], xlim=c(-1,100), main="IGRA negative")
  dev.off()

}


