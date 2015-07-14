#' Calculate a risk factor score for each patient
#'
#' \code{calcRiskFactorScore}
#'
#' @param data
#' @return data
#'

calcRiskFactorScore <- function(data){
  riskfacs <- c("WHOcut","jobrisk","numSymptoms","TBcont","PatientAge","Sex","Ethnclass","CurrHomeless","HIVpos")
  formula <- paste("TBconfirmed~", paste(riskfacs,collapse="+"), sep="")
  fit <- glm(formula, data=data, family = binomial, na.action = na.exclude)
  data$riskfacScore <- predict(fit, type = "response", na.action = na.exclude)
  data
}


#' Remove Patients Records In Data Cleaning
#'
#' \code{rmPatientsInCleaning}
#'
#' @param data
#' @return data
#'

rmPatientsInCleaning <- function(data){

  ## duplicated patients (make sure that don't lose info that we want!)
  # View(data[duplicated(data$PatientStudyID) | duplicated(data$PatientStudyID, fromLast = T),])
  data <- data[!duplicated(data$PatientStudyID),]

  ## `complete' records only before cut-off date
  data <- data[as.Date.POSIX(data$DateConsent)<=as.Date("2013-08-31") & !is.na(data$DateConsent),]

  data <- data[data$Exclude=="No" | is.na(data$Exclude) | data$Exclude=="",]

  ## non-EPTB only
  ## do we really want the cases that are _suspected_ of being PTB instead?
  data <- data[data$EPTBorPTB=="PTB" | data$EPTBorPTB=="EPTB;PTB" | is.na(data$EPTBorPTB) | data$EPTBorPTB=="",]

  ## not on our diagnostic pathway
  ### patients treated _before_ first test
  data <- data[data$preTestDrug==FALSE | is.na(data$preTestDrug),]

  ### confirmed diagnosis before first test
  data <- data[data$testDiagCon_diff>=0 | is.na(data$testDiagCon_diff),]

  data
}


#' Fill-in End Of Treatment Dates
#'
#' \code{fillInEndOfTreatmentDate}
#'
#' @param data
#' @return data
#'

fillInEndOfTreatmentDate <- function(data){
  ## add-on treatment duration time when end of treatment date not given

  TBtreatdur.num <- c("6 months"=6, "2 months"=2, "3 months"=3, "9 months"=9, "12 months"=12)
  tdur <- is.na(data$TBDrugEnd.max) & data$TBtreatdur%in%names(TBtreatdur.num)
  data$TBDrugEnd.max[tdur] <- as.Date(as.mondate(data$TBDrugStart.min)[tdur] +
                                        TBtreatdur.num[data$TBtreatdur][tdur])
  data$TBDrugEnd.max.estimate <- tdur

  data
}


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


#' prepend Not taken as extra level
#'
#' \code{joinWithLookups} takes Access database and updates one of the tables
#'
#' @param Patients single table to update the data total database
#' @return updated full database

addLevel_Nottaken <- function(data, names){
    for (i in names){
        levels(data[,i]) <- c(levels(data[,i]),"Not taken")
        data[,i] <- relevel(data[,i], ref="Not taken")
    }
    data
}


#' join data set with demographic look-up tables.
#'
#' \code{joinWithLookups} takes Access database and updates one of the tables
#'
#' @param Patients single table to update the data total database
#' @return updated full database

joinWithLookups <- function(data){

  # data(cob_lookup, envir=environment())

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


#' Assert dates are in the past
#'
#' \code{assert.dateinpast} Checks whether all of the dates all in the past
#'
#' @param data
#' @param testDate.names
#' @return NULL

assert.dateinpast <- function(data, testDate.names){
  for (i in 1:nrow(data)){
    cond <- sapply(data[i,testDate.names], function(x) difftime(x, as.Date.POSIX("2016-01-01")))
    if(!all(is.na(cond))) stopifnot(any(cond<0, na.rm=TRUE))
  }
}


#' Combine test results
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

#' Get the Date Frequencies
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



