
#' Calculate the First Test Dates for Each Patient
#'
#' Which test dates are earliest within a time window from the date of study consent.
#'
#' @param data IDEA study data set.
#' @param testDate.names Names of all test dates for each patient.
#' @param maxtime The earliest allowable date before the date of consent for a test
#'                to be taken and still be for the same episode of TB.
#' @export
#' @return testDate.min for each patient appended to original dataframe
#'
#' @examples
#'
#' testDate.names <- c("TBculttestDate","QFNtestDate", "TSPOTtestDate", "TSTtestDate", "SmeartestDate", "BALtestDate", "HistBioptestDate")
#' maxtime <- 2
#'
calc.testDate.min <- function(data,
                              testDate.names,
                              maxtime){
  require(reshape2)

  DateConsent <- data$DateConsent
  col_names <- c("PatientStudyID", testDate.names)

  res <-
    data %>%
    select(col_names) %>%
    mutate_at(vars(-PatientStudyID), as_date) %>%
    melt(id.vars = 'PatientStudyID') %>%
    group_by(variable) %>%
    mutate(DateConsent,
           keep = value - DateConsent + maxtime >= 0 & !is.na(value)) %>%
    ungroup() %>%
    filter(keep == TRUE) %>%
    group_by(PatientStudyID) %>%
    summarise(testDate.min = min(value)) %>%
    merge(data, ., by = "PatientStudyID")

  return(res)
}


#' Remove Dates That are Before the Earliest Time For the Given Suspected TB Episode
#'
#' @param data IDEA study dataset
#' @param Date.names column names to check
#' @param maxtime maimum time in window
#'
#' @inheritParams calc.Date.min
#'
#' @return original data with too early dates filled with NA
#' @export
#' @seealso \code{\link{calc.Date.min}}

rm.TooEarly <- function(data,
                        Date.names,
                        maxtime){


  whichdates <- function(x)
    (difftime(time1 = x,
              time2 = x["DateConsent"],
              units = "days") + maxtime) < 0

  tooearly <-
    data %>%
    select("DateConsent", Date.names) %>%
    apply(MARGIN = 1,
          FUN = whichdates) %>%
    t()

  data[ ,c("DateConsent", Date.names)][tooearly] <- NA

  data
}


#' Distiguish Between When Clinical Features Used to Rule-in or Rule-out TB
#'
#' Add a new column with \code{TB, other} or \code{NA}.
#'
#' @param data
#'
#' @return data

is.diagByClinicalFeatures <- function(data){

  data$clinfeatures <- rep("none", nrow(data))
  data$clinfeatures[grepl(x = data$Diagnostic_mean,
                          pattern = "Clinical") & data$TBconfirmed == "TRUE"]  <- "TB"
  data$clinfeatures[grepl(x = data$Diagnostic_mean,
                          pattern = "Clinical") & data$TBconfirmed == "FALSE"] <- "other"
  data
}


#' Calculate a risk factor score for each patient
#'
#' using a logistic regression calculate model probability predictions.
#'
#' ##TODO##
#' add clinician or hospital history covariate, maybe multilevel?
#'
#' @param data
#' @export
#' @return data with additional column

calcRiskFactorScore <- function(data){

  riskfacs <- c("WHOcut","jobrisk","numSymptoms","TBcont","PatientAge","Sex","Ethnclass","CurrHomeless","HIVpos")
  formula <- paste("TBconfirmed~", paste(riskfacs,collapse = "+"), sep = "")
  fit <- glm(formula, data = data, family = binomial, na.action = na.exclude)
  data$riskfacScore <- predict(fit, type = "response", na.action = na.exclude)
  data
}


#' Estimate Symptoms Resolved by Whether There Was a Response to Treatment or Self-Cured
#'
#' Checks if there is a response to treatment or self-cured.
#' Uses \code{all} as conservative estimates to check if the time between the observed end dates and earliest time that a patient was
#' put on treatment is smaller than the time to follow-up by clinicians.
#'
#' @param data IDEA study data set
#' @param symptomsEndDates.names
#' @param symptoms.names E.g. cough, weight loss, etc
#' @param drugReviewPeriod Time to follow-up by clinicians
#' @export
#' @return data

EstimateSymptomsResolved <- function(data,
                                     symptomsEndDates.names,
                                     symptoms.names,
                                     drugReviewPeriod){

  data$Sx_resolved <- FALSE

  for (i in 1:nrow(data)){

    cols <- symptomsEndDates.names[unlist(data[i, symptoms.names])]

    if(all(is.na(cols))){   #no symptoms
      data$Sx_resolved[i] <- FALSE
    }else{
      observedEndDates <- data[i, cols]
      observedEndDates[is.na(observedEndDates)] <- as.Date("2020/01/01")       #replace missing with large dates

      if(data$step1Diag[i]){
        data$Sx_resolved[i] <- all(sapply(observedEndDates,
                                          function(x) (x - data$TBDrugStart.min[i]) < drugReviewPeriod))      #response to treatment
        ## could alternatively have same as below?
        ## ...
      }else{
        # data$Sx_resolved[i] <- all(sapply(observedEndDates, function(x) x <= data[i,"TBculttestDate.orig"]))    #self-cured
        data$Sx_resolved[i] <- all(sapply(observedEndDates,
                                          function(x) x <= data[i, "TBculttestDate"]))
      }
    }
  }

  data
}


#' Group Together Diagnosis Outcomes
#'
#' [TODO] do over-lapping groups, not mutually exclusive
#'
#' @param data IDEA study data set
#' @export
#' @return data

GroupDiagOutcomes <- function(data){

  lookuplist <- list("LRTI"=c("LRTI",
                              "LTBI - treatment indicated",
                              "LTBI - treatment indicated;Other",
                              "LTBI - treatment indicated;Other;URTI",
                              "LTBI - treatment indicated;Other;Pneumonia",
                              "LTBI - treatment indicated;Pneumonia",
                              "LTBI - treatment indicated;URTI",
                              "LRTI;LTBI - treatment indicated;Other",
                              "LRTI;Other"),
                     "Cancer"=c("Cancer","Cancer;LTBI - treatment indicated", "Cancer;Other","Cancer;Other;Pneumonia"),
                     "Chest Infection"=c("Chest Infection", "Chest Infection;LTBI - treatment indicated", "Chest Infection;Other"),
                     "Other"=c("Other","Other;Pneumonia","Other;Sarcoidosis","Other;URTI"),
                     "Active TB"=c("Active TB", "Active TB;Other"))

  data$DiagoutcomeGrouped <- data$Diagoutcome
  lfac <- levels(data$DiagoutcomeGrouped)
  othrlevs <- lfac[!lfac %in% unlist(lookuplist)]
  x <- c(lookuplist, othrlevs)
  names(x) <- c(names(lookuplist), othrlevs)
  levels(data$DiagoutcomeGrouped) <- x

  data
}


#' Remove Patients Records In Data Cleaning
#'
#' Due to duplication, before 2013-08-31 (last date of complete records), not excluded, not EPTB only,
#' patients not treated before date of first test
#'
#' ##TODO##
#' should we remove patients whos `means of diagnosis' was IGRA?
#'
#' @param data IDEA study data set
#' @return data
#' @export

rmPatientsInCleaning <- function(data){

  ## duplicated patients (make sure that don't lose info that we want!)
  # View(data[duplicated(data$PatientStudyID) | duplicated(data$PatientStudyID, fromLast = T),])
  data <- data[!duplicated(data$PatientStudyID),]

  ## `complete' records only before cut-off date
  data <- data[as.Date.POSIX(data$DateConsent)<=as.Date("2013-08-31") & !is.na(data$DateConsent),]

  data <- subset(data, (Exclude=="No"| Exclude=="" | is.na(Exclude)))

  ## non-EPTB only
  ## do we really want the cases that are _suspected_ of being PTB instead?
  data <- subset(data, (EPTBorPTB=="PTB" | EPTBorPTB=="EPTB;PTB" | EPTBorPTB=="" | is.na(EPTBorPTB)))

  ## does it matter if we include these or not??
#   ## Heatherwood and Wrexham Park Hospital NHS Foundation Trust
#   data <- subset(data, SiteID!="H")
#
#   ## patient was TRANSFERRED TO LEICESTER FOR FU and got his study number changed to L058
#   data <- subset(data, PatientStudyID!="B048")
#
#   ## decided not to give blood, probably patient changed her mind and therefore did not ENROLL in the study
#   data <- subset(data, PatientStudyID!="B117")

  ## not on our diagnostic pathway
  ### patients treated _before_ first test
  data <- subset(data, (preTestDrug==FALSE | is.na(preTestDrug)))

  ### confirmed diagnosis before first test
  ## don't use this anymore- unreliable
  # data <- data[data$testDiagCon_diff>=0 | is.na(data$testDiagCon_diff),]


  ## withdrawn, no or invalid consent
  data <- subset(data, !PatientStudyID%in%c("B008",  "N180",  "N301",  "G045",  "G057",  "G059",  "C026",  "C032",   "C033",  "F017",  "F047",  "F050",  "F057"))


  data
}


#' Fill-in End Of Treatment Dates
#'
#' If we have the duration of treatment and start of treatment but no
#' end date then we can calculate the end of treatment date.
#'
#' @param data IDEA study data set
#' @return data
#' @export

fillInEndOfTreatmentDate <- function(data){
  ## add-on treatment duration time when end of treatment date not given

  TBtreatdur.num <- c("6 months"=6, "2 months"=2, "3 months"=3, "9 months"=9, "12 months"=12)
  tdur <- is.na(data$TBDrugEnd.max) & data$TBtreatdur%in%names(TBtreatdur.num)
  data$TBDrugEnd.max[tdur] <- as.Date(as.mondate(data$TBDrugStart.min)[tdur] +
                                        TBtreatdur.num[data$TBtreatdur][tdur])
  data$TBDrugEnd.max.estimate <- tdur

  data
}


#' Read From Database Extract
#'
#' \code{joinLevels} takes Access database and updates one of the tables.
#'
#' @param Patients Single table to update the data total database
#' @return codes Updated full database
#' @export

joinLevels <- function (codes, lookuplist){
    if(!is.factor(codes)){break}
    lfac <- levels(codes)
    othrlevs <- lfac[!lfac %in% unlist(lookuplist)]
    x <- c(lookuplist, othrlevs)
    names(x) <- c("UK", othrlevs)
    levels(codes) <- x
    codes
}


#' Prepend Not Taken as Extra Level
#'
#' \code{joinWithLookups} takes Access database and updates one of the tables.
#'
#' @param Patients single table to update the data total database
#' @return data Updated full database
#' @export
addLevel_Nottaken <- function(data, names){
    for (i in names){
        levels(data[,i]) <- c(levels(data[,i]),"Not taken")
        data[,i] <- relevel(data[,i], ref="Not taken")
    }
    data
}


#' Join Data Set With Demographic Look-up Tables
#'
#' \code{joinWithLookups} takes Access database and updates one of the tables.
#'
#' @param Patients single table to update the data total database
#' @return data Updated full database
#' @export
joinWithLookups <- function(data){

  # data(cob_lookup, envir=environment())

  cob_lookup$Region <- toupper(cob_lookup$Region)

  data <- merge(data, ethn_lookup[,c("EthnID", "Ethnclass", "Ethndetail")], by.x="Ethn", by.y="EthnID", all.x=TRUE)
  data <- merge(data, site_lookup, all.x=TRUE)
  data <- merge(data, cob_lookup, by.x="CtryBirth", by.y="CountryrecordID", all.x=TRUE)

  data$WHOcut <- cut(data$TBincidence, breaks=c(0,40,100,150,200,400,10000), right=FALSE)

  data
}


#' Convert Different Classes to Date Class
#'
#' \code{as.Date.POSIX} takes Access database and updates one of the tables.
#'
#' @param Patients single table to update the data total database
#' @return data Updated full database
#' @export
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
#' @export
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
#' @export
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

