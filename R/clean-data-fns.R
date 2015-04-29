#' Read from database extract.
#'
#' \code{joinLevels} takes Access database and updates one of the tables
#'
#' @param Patients single table to update the data total database
#' @return updated full database
#'

joinLevels <- function (codes, lookuplist)
{
    if(!is.factor(codes)){break}
    lfac <- levels(codes)
    othrlevs <- lfac[!lfac %in% unlist(lookuplist)]
    x <- c(lookuplist, othrlevs)
    names(x) <- c("UK", othrlevs)
    levels(codes) <- x
    codes
}


#' join data set with demographic look-up tables.
#'
#' \code{joinWithLookups} takes Access database and updates one of the tables
#'
#' @param Patients single table to update the data total database
#' @return updated full database

joinWithLookups <- function(data){
  #

  cob_lookup$Region <- toupper(cob_lookup$Region)

  data <- merge(data, ethn_lookup[,c("EthnID", "Ethnclass", "Ethndetail")], by.x="Ethn", by.y="EthnID", all.x=TRUE)
  data <- merge(data, site_lookup, all.x=TRUE)
  data <- merge(data, cob_lookup, by.x="CtryBirth", by.y="CountryrecordID", all.x=TRUE)

  data$WHOcut <- cut(data$TBincidence, breaks=c(0,40,100,150,200,400,10000), right=FALSE)

  data
}


#' Convert different classes to Date class.
#'
#' \code{as.Date.POSIX} takes Access database and updates one of the tables
#'
#' @param Patients single table to update the data total database
#' @return updated full database

as.Date.POSIX <- function(date){
  # http://stackoverflow.com/questions/8788817/r-converting-posixct-dates-with-bst-gmt-tags-using-as-date
    ## coerce to class Date

    if(is.character(date)){
        date <- as.Date(date)
    }else if(is.factor(date)){
        date <- as.Date(as.character(date), "%d/%m/%Y")
    }else if(is.POSIXt(date)){
        date <- as.Date(date+3600)
    }

    date
}


#' Combine test results.
#'
#' \code{combineTestResults} returns joined tables of inputs and predicted outcomes
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories
#'
combineTestResults <- function(test1, test2){

    test1 <- as.character(test1)
    test2 <- as.character(test2)
    test1[is.na(test1)] <- ""
    test2[is.na(test2)] <- ""
    res <- NULL

    for(i in 1:length(test1)){

        ## both not empty
        if(!test1[i]=="" & !test2[i]==""){

            if(test1[i]=="NEGATIVE" & test2[i]=="NEGATIVE"){
                res[i] <- "NEGATIVE"}
            else if(test1[i]=="POSITIVE" | test2[i]=="POSITIVE"){
                res[i] <- "POSITIVE"}
        }else{
            res[i] <- paste(test1[i], test2[i], sep="")}
    }
    res[res==""] <- NA

    res
}

#' getDateFrequencies
#'
#' \code{getDateFrequencies} returns joined tables of inputs and predicted outcomes
#'
#' @param data Individual patient records
#' @param allgrids The idealised pathway look-up tables
#' @return The joined tables with aggregated outcome wrt Dosanjh categories

getDateFrequencies <- function(colnames, data){
    ## unique dates and frequencies on those dates
    ## e.g. for starting treatment

    df <- matrix(data=NA, nrow=nrow(data), ncol=length(colnames)*2)
    colnames(df) <- c(paste("DateUnique", 1:length(colnames), sep=""), paste("Freq", 1:length(colnames), sep=""))

    for (i in 1:nrow(df)){
        tab <- table(t(data[i,colnames]))
        if (length(tab)!=0){
            df[i,1:length(tab)] <- names(tab)
            df[i,length(colnames)+1:length(tab)] <- as.numeric(tab)
        }
    }
df
}



