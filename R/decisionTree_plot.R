#' Plot data paths of decision tree with frequencies
#'
#' \code{decisionTree.plot}
#'
#' @param data casted frequency along each pathway
#' @param ysteps which branching pattern
#' @return none

decisionTree.plot <- function(data, ysteps="exp", out=NA, newplot=TRUE, col="grey"){
  ## data: casted frequeny along each pathway
  ## ystep: which branching pattern

  require(plyr)

  if(!is.na(out)) colnames(data)[colnames(data)==out] <- "out"

  data.branches <- data[,-ncol(data)]

  entries <- c(T,F)
  all <- expand.grid(data.frame(matrix(entries, nrow=2, ncol=ncol(data.branches))))
  names(all) <- names(data.branches)
  data <- merge(all,data, all.x = TRUE)

  ## six branches
  scale <- 50
  shift <- 0.8
  pow <- 0.85

  x <- 0:ncol(data.branches)
  # seq_y <- c(0, shift*scale, 1/(1:(length(x)-2))^2*scale)
  seq_y <- c(0, shift*scale, exp(-1:-(length(x)-2))^pow*scale)

  posneg <- c(1,-1)
  altsign <- data.frame(0,expand.grid(data.frame(matrix(posneg, nrow=2, ncol=ncol(data.branches)))))

  y <- t(apply(altsign, 1, function(x) cumsum(seq_y*x)))

  if(newplot) plot(1,1, type="n", xlim=c(min(x)-0.5,max(x)+0.5), ylim=c(min(y)-2,max(y)+5), axes=FALSE, xlab="", ylab="")

  s <- seq(length(x)-1)

  rulelist <- lapply(apply(data[,-ncol(data)],2,list),unlist)

  linewidth <- list()
  linewidth[[1]] <- aggregate(data$out, by=list(data[,1]), FUN=function(x) sum(x, na.rm=T)*(1+runif(1,max=0.001)))

  for (i in 2:(ncol(data)-1)){
    linewidth[[i]] <- aggregate(data$out, by=rulelist[1:i], FUN=function(x) sum(x, na.rm=T)*(1+runif(1,max=0.001)))
    linewidth[[i]][,ncol(linewidth[[i]])][linewidth[[i]][,ncol(linewidth[[i]])]==0] <- NA
  }

  y_stacked <- list()
  for (i in 1:length(linewidth)){
    y_stacked[[i]] <- data.frame(xcoord=x[i], y, linewidth[[i]])
  }

  y_stacked <- ldply(y_stacked, data.frame)
  y_stacked$x[y_stacked$x==0] <- NA
  # y_stacked <- y_stacked[!(duplicated(y_stacked$x) & y_stacked$x>15 | is.na(y_stacked$x)),]
  y_stacked <- y_stacked[!duplicated(y_stacked$x),]

  for (i in 1:nrow(y_stacked)){

    ## diagonal joining lines
#     segments(y_stacked$xcoord[i],
#              y_stacked[i, y_stacked$xcoord[i]+2],
#              y_stacked$xcoord[i]+1,
#              y_stacked[i, y_stacked$xcoord[i]+3],
#              col=col,
#              lwd=y_stacked$x[i])

    ## right-angled joining lines
    segments(y_stacked$xcoord[i],
             y_stacked[i, y_stacked$xcoord[i]+2],
             y_stacked$xcoord[i],
             y_stacked[i, y_stacked$xcoord[i]+3],
             col=col,
             lwd=y_stacked$x[i], lend="square")
    segments(y_stacked$xcoord[i],
             y_stacked[i, y_stacked$xcoord[i]+3],
             y_stacked$xcoord[i]+1,
             y_stacked[i, y_stacked$xcoord[i]+3],
             col=col,
             lwd=y_stacked$x[i], lend="butt")
  }

  for (i in 1:nrow(y)){
    points(x,y[i,], col="lightblue", pch=19, cex=1)
  }

  # add text and labels
  text(x[-1]-0.5,max(y)+5,names(data.branches))


}


#' Plot multiple data paths of decision tree with frequencies
#'
#' \code{decisionTreeMultiple.plot}
#'
#' @param dt casted frequency along each pathway
#' @param names which branching pattern
#' @param outcol
#' @param ncols to use 4 columns or 2
#' @return none

decisionTreeMultiple.plot <- function(dt, names, outcol, ncols=4){

  if(!is.vector(names)) stop("names is not character vector")

  trsp <- 0.4
  cols <- c(rgb(0,1,0,trsp), rgb(0,0,1,trsp), rgb(1,0,1,trsp), rgb(1,0,0,trsp))

  ## split by Dosangh category
  ### overlayed
  par(mfrow=c(1,1))
  data <- data.frame(dt[,names], logw=2*log(as.numeric(dt1[,outcol])+1))
  decisionTree.plot(data, out="logw", col=cols[1])

  for (i in 1:(ncols-1)){
    data <- data.frame(dt[,names], logw=2*log(as.numeric(dt1[,outcol+i])+1))
    decisionTree.plot(data, out="logw", newplot=F, col=cols[i+1])
  }


  ### gridded
  par(mfrow=c(ceiling(ncols/2),2))
  par(mar=c(0,0,0,0))
  par(oma=c(0,0,3,0))

  data <- data.frame(dt[,names], logw=log(as.numeric(dt1[,outcol])+1))
  decisionTree.plot(data, out="logw")

  for (i in 1:(ncols-1)){
    data <- data.frame(dt[,names], logw=log(as.numeric(dt1[,outcol+i])+1))
    decisionTree.plot(data, out="logw")
  }
}



