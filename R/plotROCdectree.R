
plotROCdectree <- function(tab){

  ## dark blue: don't treat
  ## light blue: treat
  ## green: ideal total
  ## red: data total
  ## black: smear positive

  transf <- function(x) log(as.numeric(x))*2 #x/20

  trsp <- 0.3
  cols <- c(rgb(0,1,0,trsp), rgb(0,0,1,trsp), rgb(1,0,1,trsp), rgb(1,0,0,trsp))


  ## observations ##

  obs <- tab[[1]]$treat1obs
  nobs.original <- nrow(obs)

  addTestResults <- function(obs){

    obs.naomit <- na.omit(obs)

    ntests <- 4       # hard coded 4 tests
    sumTestResults <- function(testname){
      cbind(names(table(as.character(obs.naomit[,testname]))),
            sapply((ntests+1):ncol(obs), function(x)
              aggregate( as.numeric(obs.naomit[,x]), by=list(obs.naomit[,testname]), FUN=sum)$x))
    }


    ## aggragate to sum over one test results only

    x <- sapply(1:ntests, sumTestResults)
    sumTestRows <- NULL
    for (i in 1:ntests){
      temp <- matrix(data=NA, nrow = nrow(x[[i]]), ncol=ntests)
      temp[,i] <- x[[i]][,1]
      sumTestRows <- rbind(sumTestRows, cbind(temp, x[[i]][,-1]))
    }
    colnames(sumTestRows) <- names(obs)


    ## idealised treatment decision with one test only

    DosanjhSums <- t(apply(sumTestRows, 1, function(x) as.numeric(x[1:4+ntests]) + as.numeric(x[1:4+(2*ntests)])))
    DosanjhSums <- cbind(melt(sumTestRows[,1:4], measure.vars = colnames(sumTestRows)[1:4], na.rm=T), DosanjhSums)
    zeros <- rep(0,ntests)
    DsumwithZeros <- matrix(NA, nrow=nrow(DosanjhSums), ncol=2*ntests)
    i <- 1
    for (value in DosanjhSums$value){

      if(value%in%c("NEGATIVE","Not taken")) DsumwithZeros[i,] <- c(zeros,as.numeric(DosanjhSums[i,1:ntests+3]))
      if(value=="POSITIVE")  DsumwithZeros[i,] <- c(as.numeric(DosanjhSums[i,1:ntests+3]), zeros)
      if(value=="TRUE") DsumwithZeros[i,] <- c(as.numeric(DosanjhSums[i,1:ntests+3]), zeros)
      if(value=="FALSE")  DsumwithZeros[i,] <- c(zeros,as.numeric(DosanjhSums[i,1:ntests+3]))
      i <- i+1
    }

    singleTestFreq <- sumTestRows
    singleTestFreq[,(ntests+1):(3*ntests)] <- DsumwithZeros

    rbind(obs, sumTestRows, singleTestFreq)
  }

  obs <- addTestResults(obs)

  obs <- transform(obs,
                   "X12_TRUE"=as.numeric(X1_TRUE)+as.numeric(X2_TRUE),
                   "X12_FALSE"=as.numeric(X1_FALSE)+as.numeric(X2_FALSE))
  obs <- transform(obs,
                   TPrate=as.numeric(X12_TRUE)/(as.numeric(X12_TRUE)+as.numeric(X12_FALSE)),
                   FPrate=as.numeric(X4_TRUE)/(as.numeric(X4_TRUE)+as.numeric(X4_FALSE)))
  obs <- transform(obs,
                   TPsd=sqrt(TPrate*(1-TPrate)/(as.numeric(X12_TRUE)+as.numeric(X12_FALSE))),
                   FPsd=sqrt(FPrate*(1-FPrate)/(as.numeric(X4_TRUE)+as.numeric(X4_FALSE))))
  # x$FPrate[is.na(x$FPrate)] <- 0
  # x$TPrate[is.na(x$TPrate)] <- 0


  ## idealised ##

  ideal <- tab[[1]]$treat1est

  ideal <- addTestResults(ideal)

  ideal <- transform(ideal,
                     "X12_TRUE"=as.numeric(X1_TRUE)+as.numeric(X2_TRUE),
                     "X12_FALSE"=as.numeric(X1_FALSE)+as.numeric(X2_FALSE))
  ideal <- transform(ideal,
                     TPrate=as.numeric(X12_TRUE)/(as.numeric(X12_TRUE)+as.numeric(X12_FALSE)),
                     FPrate=as.numeric(X4_TRUE)/(as.numeric(X4_TRUE)+as.numeric(X4_FALSE)))
  ideal <- transform(ideal,
                     TPsd=sqrt(TPrate*(1-TPrate)/(as.numeric(X12_TRUE)+as.numeric(X12_FALSE))),
                     FPsd=sqrt(FPrate*(1-FPrate)/(as.numeric(X4_TRUE)+as.numeric(X4_FALSE))))


  ideal.treatTF <- rowSums(ideal[,c("TPrate", "FPrate")],na.rm=T)>0


  ## plots ##

  par(mfrow=c(1,1))

  plot(obs$FPrate[2:nobs.original-1], obs$TPrate[2:nobs.original-1],
       cex=transf(obs$Total)[2:nobs.original-1],
       xlim=c(0,0.2), pch=19, col=cols[ideal.treatTF[2:nobs.original-1]+1],
       xlab="Proportion of non-TB presumptively treated", ylab="Proportion TB presumptively treated")
  # xlab="False positive rate", ylab="True positive rate")
  points(obs$FPrate[1], obs$TPrate[1], cex=transf(obs$Total[1]), lwd=2, pch=19, col="red")
  points(obs$FPrate[nobs.original], obs$TPrate[nobs.original], lwd=2, cex=transf(obs$Total[nobs.original]))
  abline(a=0, b = 1, lty=2)
  text(obs$FPrate-0.001, obs$TPrate+0.01, 1:nrow(obs), cex=0.7)


  # points(x$FPrate[2:nrow(x)-1]+rnorm(n=nrow(x)-1,sd = 0.005), x$TPrate[2:nrow(x)-1]+rnorm(n=nrow(x)-1,sd = 0.005), cex=transf(x$Total[2:nrow(x)-1]),
  # pch=3)#, col=cols[1])#, type = n)
  points(ideal$FPrate[nobs.original], ideal$TPrate[nobs.original], lwd=2, cex=transf(ideal$Total[nobs.original]), pch=3)

  points(obs$FPrate[nobs.original:nrow(obs)], obs$TPrate[nobs.original:nrow(obs)],
         lwd=2, cex=transf(obs$Total[nobs.original:nrow(obs)]), pch=2)


  ## add ellipse or vertical and horizontal lines
  ## indicating variation of FP rate and TP rate
  ## e.g. \sqrt(p*(1-p)/n))
  ## http://stackoverflow.com/questions/9231702/ggplot2-adding-two-errorbars-to-each-point-in-scatterplot
  ##http://stats.stackexchange.com/questions/29641/standard-error-for-the-mean-of-a-sample-of-binomial-random-variables


  df <- data.frame(x = c(obs$FPrate, ideal$FPrate[nobs.original]),
                   y = c(obs$TPrate, ideal$TPrate[nobs.original]),
                   ymin = c(obs$TPrate - obs$TPsd, ideal$TPrate[nobs.original]-ideal$TPsd[nobs.original]),
                   ymax = c(obs$TPrate + obs$TPsd, ideal$TPrate[nobs.original])+ideal$TPsd[nobs.original],
                   xmin = c(obs$FPrate - obs$FPsd, ideal$FPrate[nobs.original]-ideal$FPsd[nobs.original]),
                   xmax = c(obs$FPrate + obs$FPsd, ideal$FPrate[nobs.original]+ideal$FPsd[nobs.original]))


#   colplot <- ideal.treatTF+4
#   colplot[1] <- 1
#   colplot[length(colplot)] <- 2
#   colplot <- c(colplot,3)
  colplot <- ideal.treatTF+4
  colplot[1] <- 1
  colplot[nobs.original:length(colplot)] <- 0
  colplot <- c(colplot,3)

  fig <- ggplot(data = df, aes(x = x,y = y)) + theme_bw() + theme(legend.position="none") +
    xlim(0,0.3) +
    geom_point(aes(size=4, colour=colplot)) +
    geom_errorbar(aes(ymin = ymin,ymax = ymax, colour=colplot, size=2, alpha=0.3)) +
    geom_errorbarh(aes(xmin = xmin,xmax = xmax, colour=colplot, size=2, alpha=0.3)) +
    xlab("Proportion of non-TB presumptively treated") +
    ylab("Proportion TB presumptively treated") +
    geom_text(label=1:length(colplot), hjust=1, vjust=-1) +
    scale_colour_gradientn(colours=rainbow(5))
  print(fig)


  df <- data.frame(x = c(obs$X12_FALSE, ideal$X12_FALSE[nrow(ideal)]),
                   y = c(obs$X12_TRUE, ideal$X12_TRUE[nrow(ideal)]))
  df2 <- data.frame(x = c(as.numeric(obs$X4_FALSE), as.numeric(ideal$X4_FALSE[nrow(ideal)])),
                    y = c(as.numeric(obs$X4_TRUE), as.numeric(ideal$X4_TRUE[nrow(ideal)])))

  par(mfrow=c(1,2))
  plot(df$x, df$y,
       xlab="Number of non-TB cases presumptively treated",
       ylab="Number of TB cases presumptively treated",
       col=rainbow(5)[colplot], pch=19)#, log="xy")
  text(df$x-1, df$y+1, 1:nrow(df), cex=0.7)
  plot(df2$x, df2$y,
       xlab="Number of TB cases not presumptively treated",
       ylab="Number of non-TB cases not presumptively treated",
       col=rainbow(5)[colplot], pch=19)#, log="xy")
  text(df2$x-0.5, df2$y+0.5, 1:nrow(df2), cex=0.7)

}











