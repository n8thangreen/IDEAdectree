

## bootstrap sampled patients
sample.IDs <- function(cat){
  sapply(numNotRuledOut.new[[cat]], function(x)
    sample(data$PatientStudyID[data$DosanjhGrouped==cat], size=x, replace = FALSE), simplify = F)
}


calc.ruleOutDiag <- function(cat, N, FNtime, ruleouttime){
  function(nruledin){
    ifelse(cat==4,
           return(ruleouttime*N),
           return((ruleouttime*N) + ((N-nruledin)*FNtime))
    )}
}

calc.daysToDiag <- function(cat){

  N <- NumDosanjh[cat]
  ruleOutDiag <- calc.ruleOutDiag(cat, N, FNtime=0, ruleouttime=1)

  sapply(sampledIDs[[cat]], function(x){
    ((mean(data$start.to.diag[data$PatientStudyID%in%x], na.rm=T)*length(x)) + ruleOutDiag(length(x)))/N
  })
}

calc.ruleOutCost <- function(cat, N, FNcost){
  function(ruleoutcost, nruledin){
    ifelse(cat==4,
           return(ruleoutcost*N),
           return((ruleoutcost*N) + ((N-nruledin)*FNcost))
    )}
}

calc.meancost.new <- function(cat){

  ruleOutCost <- calc.ruleOutCost(cat, N=NumDosanjh[cat], FNcost=0)
  mapply(function(x,y) ((mean(data$totalcost[data$PatientStudyID%in%x], na.rm=T)*length(x)) + ruleOutCost(y, length(x)))/NumDosanjh[cat],
         sampledIDs[[cat]], testcost)
}


##############
## plotting ##
##############
# http://stackoverflow.com/questions/5044678/stack-contour-plot-in-r  #overlaying contours

plot.surface_smooth <- function(out){

  axis <- 1:(length(out[[1]]$Cruleout_hat)/length(unique(out[[1]]$testcost)))  # unique sens/spec
  d <- data.frame(Sensitivity=out[[1]]$sensitivity[axis], Specificity=out[[1]]$specificity[axis], C=out[[1]]$Cruleout_hat[axis])
  x <- matrix(out[[1]]$Cruleout_hat[axis], ncol=length(unique(out[[1]]$sensitivity[axis])), nrow=length(unique(out[[1]]$sensitivity[axis])))

  ## effective combined, single test
  # d <- data.frame(Sensitivity=out[[1]]$sensitivity.effective[axis], Specificity=out[[1]]$specificity.effective[axis], C=out[[1]]$Cruleout_hat[axis])
  d <- d[!duplicated(d),]

  ## lattice pkg
  #   image(x=rev(thresh), y=rev(thresh), z=x)
  #   contour(x=rev(thresh), y=rev(thresh), z=x[], add = T)
  # lattice::levelplot(x)
  # lattice::contourplot(x,)

  # ggplot(d, aes(sens, spec, fill=C)) + geom_tile()+scale_fill_gradient(limits = c(100, 500), low = "yellow", high = "red")
  ggplot(d, aes(Sensitivity,Specificity,z=C)) + geom_tile(aes(fill=C))+
    # stat_contour(bins=6, aes(Sensitivity,Specificity,z=C), color="black", size=0.6)+
    stat_contour(aes(Sensitivity,Specificity,z=C, colour="..level.."), color="black", size=0.6, breaks=seq(-1000,5000,by=20))+
                   # c(seq(-100,-0,by=50),seq(-1,1,by=50), seq(0,1000,by=50)))+
    scale_fill_gradientn(colours=brewer.pal(6,"YlOrRd")) + theme_bw() + ylim(0.5,1)#, limits=c(40,100)
  # direct.label(v)
}

plot.surface_solid <- function(out){
  # http://stackoverflow.com/questions/10981324/ggplot2-heatmap-with-colors-for-ranged-values

  axis <- 1:(length(out[[1]]$Cruleout_hat)/length(unique(out[[1]]$testcost)))  # unique sens/spec
  mat <- data.frame(Sensitivity=out[[1]]$sensitivity[axis], Specificity=out[[1]]$specificity[axis], C=out[[1]]$Cruleout_hat[axis])
  len <- 8
  brks <- cut(mat$C,breaks=seq(0,200,len=len))
  brks <- gsub(","," - ",brks,fixed=TRUE)
  mat$brks <- gsub("\\(|\\]","",brks)  # reformat guide labels

  ggplot(mat,aes(Sensitivity, Specificity))+
    geom_tile(aes(fill=brks))+
    scale_fill_manual("Z", values=brewer.pal(len,"YlOrRd"))+
    # scale_fill_manual(values = colorRampPalette(c("orange", "yellow"))(len))+   #if len>9
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    coord_fixed()
}


plot.surface_threshold <- function(out1, out2, threshold=90){

  axis <- 1:(length(out1[[1]]$Cruleout_hat)/length(unique(out1[[1]]$testcost)))  # unique sens/spec
  d1 <- data.frame(Sensitivity=out1[[1]]$sensitivity[axis], Specificity=out1[[1]]$specificity[axis], C=out1[[1]]$Cruleout_hat[axis])
  d2 <- data.frame(Sensitivity=out2[[1]]$sensitivity[axis], Specificity=out2[[1]]$specificity[axis], C=out2[[1]]$Cruleout_hat[axis])
  d1$C<-as.numeric(d1$C>threshold)
  d1$C <- d1$C + as.numeric(d2$C>threshold)

  ggplot(d1, aes(Sensitivity,Specificity,z=C)) + geom_tile(aes(fill=C)) +
    scale_fill_gradientn(colours=brewer.pal(6,"Greys")) + theme_bw() #, limits=c(40,100)
}



