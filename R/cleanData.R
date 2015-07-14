#' Clean IDEA clinical data.
#'
#' \code{cleanData} is a high-level function to clean IDEA clinical data for analysis
#'
#' @param data Individual patient records
#' @return data

cleanData <- function(data){
  ## data: file address (string)
  ##

  # load(file=data)


  # TB tests ------------------------------------------------------------------

  names(data) <- gsub("testRes", "", names(data))

  # Discretised Groups -----------------------------------------------
  data$TSTcut <- cut(data$TST, breaks=c(0,6,15,100), right=FALSE)  #used by Yemesi


  test.names <- c("QFN", "TSPOT", "TST", "TSTcut", "Smear", "TBcult", "CSF", "BAL",
                  "HistBiop", "NeedleAsp", "PCR", "CXR", "CT", "MRI")

  data[,c(test.names,"Othertest1Res","Othertest2Res")] <- data.frame(apply(data[,c(test.names,"Othertest1Res","Othertest2Res")], 2, toupper))

  data$IGRA <- combineTestResults(data$QFN, data$TSPOT)

  data$imaging <- combineTestResults(data$CT, data$MRI)


  # Combine with BCG status -------------------------------------------------------------

  data$TSTres <- ""
  data$TSTres[!is.na(data$TST) & data$TSTcut!="[0,6)" & data$PrevBCG==FALSE] <- "POSITIVE"
  data$TSTres[!is.na(data$TST) & data$TSTcut=="[15,100)" & data$PrevBCG==TRUE] <- "POSITIVE"
  data$TSTres[!is.na(data$TST) & data$TSTcut=="[0,6)"] <- "NEGATIVE"
  data$TSTres[!is.na(data$TST) & data$TSTcut=="[6,15)" & data$PrevBCG==TRUE] <- "NEGATIVE"

  data$IGRA[data$IGRA=="BORDERLINE"] <- NA
  data$IGRA[data$IGRA=="INDETERMINATE"] <- NA

  data$IGRAorTST <- combineTestResults(data$TSTres, data$IGRA)


  # look for PET in Othertests---------------------------------------------

  PETtest1taken <- grepl("PET", data$Othertest1Name, ignore.case=TRUE)
  PETtest2taken <- grepl("PET", data$Othertest2Name, ignore.case=TRUE)

  ## PET test results, when taken
  PET1 <- PET2 <- rep(NA, nrow(data))
  PET1[PETtest1taken] <- as.character(data$Othertest1Res[PETtest1taken])
  PET2[PETtest2taken] <- as.character(data$Othertest2Res[PETtest2taken])

  data$PET <- combineTestResults(PET1, PET2)

  ## PET test date, when taken. When both dates use test 2 (?)
  PETtestDate <- rep(as.Date("1900-01-01"), nrow(data))
  PETtestDate[PETtest1taken] <- as.Date.POSIX(data$Othertest1Date[PETtest1taken])
  PETtestDate[PETtest2taken] <- as.Date.POSIX(data$Othertest2Date[PETtest2taken])

  data$PETtestDate <- PETtestDate

  test.names <- c("IGRA", "imaging", test.names, "PET", "TSTres", "IGRAorTST")



  data <- joinWithLookups(data)


  # Dates -------------------------------------------------------------------

  TBdrugStart.names <- c("TBdrug1Start", "TBdrug2Start", "TBdrug3Start", "TBdrug4Start", "TBdrug5Start", "TBdrug6Start")
  TBDrugEnd.names <- c("TBDrug1End", "TBDrug2End", "TBDrug3End", "TBDrug4End", "TBDrug5End", "TBDrug6End")
  testDate.names <- c("TBculttestDate","QFNtestDate", "TSPOTtestDate", "TSTtestDate", "SmeartestDate", "BALtestDate", "HistBioptestDate",
                      "NeedleAsptestDate", "PCRtestDate", "CXRtestDate", "CTtestDate", "MRItestDate", "Othertest1Date", "Othertest2Date","PETtestDate")

  for (i in c(TBdrugStart.names, TBDrugEnd.names, testDate.names,
              "EntryUK", "DateDiagCon", "DateConsent")){
    data[,i] <- as.Date.POSIX(data[,i])
  }

  ## missing dates
  data[data=="1900-01-01"] <- NA

  data$EntryUK_year <- year(data$EntryUK)

  ## how can we know which are filled-in dates and which are the accurate ones?
  ## take the hit of omitting info for true dates too
  ## not so bad for day-month but for day this is not good
  ### ie all Junes or all 1st of month dates
  for (i in 1:ncol(data)){
    data[grepl("-06-15", data[,i]),i] <- NA}


  assert.dateinpast(data, testDate.names)


  ## combine drug start dates for multiple patient records
  for (id in unique(data$PatientStudyID)[!is.na(unique(data$PatientStudyID))]){

    for (i in TBdrugStart.names){
      data[data$PatientStudyID==id & !is.na(data$PatientStudyID), i] <-
        min(data[data$PatientStudyID==id & !is.na(data$PatientStudyID), i], na.rm=T)
    }
    for (j in TBDrugEnd.names){
      data[data$PatientStudyID==id & !is.na(data$PatientStudyID), j] <-
        min(data[data$PatientStudyID==id & !is.na(data$PatientStudyID), j], na.rm=T)
    }
  }


  data$TBDrugStart.min <- apply(data[,TBdrugStart.names], 1, min, na.rm=T)
  data$TBDrugEnd.max   <- apply(data[,TBDrugEnd.names], 1, max, na.rm=T)
  data$testDate.min    <- apply(data[,testDate.names[-1]], 1, min, na.rm=T)   #exclude culture test date in this

  data$TBDrugStart.min <- as.Date.POSIX(data$TBDrugStart.min)
  data$TBDrugEnd.max   <- as.Date.POSIX(data$TBDrugEnd.max)
  data$testDate.min    <- as.Date.POSIX(data$testDate.min)

  data$DateVisitFU <- as.Date.POSIX(data$DateVisitFU)


  # lower limits culture report dates -----------------------------------------

  dur <- c(POSITIVE=1, NEGATIVE=42)   #days i.e. 6 weeks
  data$TBculttestDate.orig <- data$TBculttestDate
  dur.each <- as.vector(dur[as.character(data$TBcult)])
  data <- transform(data, TBculttestDate.resMin = testDate.min+dur.each)
  data$TBculttestDate <- apply(data, 1,
                               function(x) max(x["TBculttestDate"], x["TBculttestDate.resMin"], na.rm=T))
  data$TBcultCens <- data$TBculttestDate!=data$TBculttestDate.orig


  # time-to-events ----------------------------------------------------------

  calcTimeToEvent <- function(testDate) difftime(testDate, data$testDate.min, units="days")


  data <- transform(data,
                    TBconfirmed = (Diagoutcome%in%c("Active TB", "Active TB;Other")),

                    preCultDrug = (TBculttestDate>TBDrugStart.min),
                    preTestDrug = (testDate.min>TBDrugStart.min),

                    DrugCult_diff = round(difftime(TBculttestDate, TBDrugStart.min, units="days")),
                    testDrug_diff = calcTimeToEvent(TBDrugStart.min),
                    testDiagCon_diff = round(calcTimeToEvent(DateDiagCon)),
                    testCult_diff = calcTimeToEvent(TBculttestDate),
                    testCultorig_diff = calcTimeToEvent(TBculttestDate.orig),
                    EntryUKtest_diff = year(testDate.min) - EntryUK_year,

                    start.to.Smear = calcTimeToEvent(SmeartestDate),
                    start.to.HistBiop = calcTimeToEvent(HistBioptestDate),
                    start.to.TBcultorig = calcTimeToEvent(TBculttestDate.orig),
                    start.to.BAL = calcTimeToEvent(BALtestDate),
                    start.to.PCR = calcTimeToEvent(PCRtestDate),
                    start.to.TST = calcTimeToEvent(TSTtestDate),
                    start.to.CT = calcTimeToEvent(CTtestDate),
                    start.to.MRI = calcTimeToEvent(MRItestDate),
                    start.to.PET = calcTimeToEvent(PETtestDate),
                    start.to.QFN = calcTimeToEvent(QFNtestDate),
                    start.to.TSPOT = calcTimeToEvent(TSPOTtestDate),

                    start.to.FU = calcTimeToEvent(DateVisitFU)
  )

  # View(data.frame(data$PatientStudyID, data$testDate.min, data$DateVisitFU, data$DateVisitFU0, data$start.to.FU)[order(data$PatientStudyID),])  #check

  data$start.to.Imaging <- pmax(data$start.to.CT, data$start.to.MRI, data$start.to.PET, na.rm=T)
  data$start.to.IGRA <- pmax(data$start.to.QFN, data$start.to.TSPOT, na.rm=T)

  ## only interested in 2 month followup
  # data$start.to.FU[data$VisitFU!="2 month FU"] <- NA

  data$preTestDrug[is.na(data$preTestDrug)] <- FALSE

  data$step1Diag <- (data$preCultDrug & !data$preTestDrug)
  data$step1Diag[is.na(data$step1Diag)] <- FALSE


  ## remove previous TB and LTBI patients
  ##TODO##
  ## whats the field for this??


  data <- fillInEndOfTreatmentDate(data)

  data$TBDrug_diff <- difftime(data$TBDrugEnd.max, data$TBDrugStart.min, units="days")

  ## this isn't perfect because there may be treatment gaps but is an ok approximation
  drugReviewPeriod <- 63     #i.e. ~2 months
  data$treatResponse <- (as.numeric(data$TBDrug_diff) > drugReviewPeriod)
  data$testDrug_diff_plus63days <- data$testDrug_diff + drugReviewPeriod


  TBdrugStart.freq <- getDateFrequencies(TBdrugStart.names, data)
  TBDrugEnd.freq <- getDateFrequencies(TBDrugEnd.names, data)
  testDate.freq <- getDateFrequencies(testDate.names, data)


  data <- estimateTimeToDiagnosis(data)


  # BAL ---------------------------------------------------------------------

  ## split Smear and Culture in to 2 types BAL/non-BAL
  data$SmearBAL <- data$Smear%in%c("BAL","bal")
  data$SmearNonBAL <- !data$SmearBAL
  data$TBcultBAL <- data$TBcultsite%in%c("BAL","bal")
  data$TBcultnonBAL <- !data$TBcultBAL


  ##TODO##
  ## what do we do with the BAL column, because we don't know if its culture or smear?

  # data$SmearNonBAL[data$TypeSmear%in%c("BAL","bal")] <- NA
  # data$SmearBAL[is.na(data$BAL) & data$TypeSmear%in%c("BAL","bal")] <- data$TypeSmear[is.na(data$BAL) & data$TypeSmear%in%c("BAL","bal")]

  ## also BAL could be in the Othertest fields...


  ##TODO##
  ## only include tests that report _before_ culture reports
  ## if they're after then substitute in "Not taken"
  # names(data)[grep("testDate", names(data))]


  ## relevel factors -----------------------------------------------------------------------
  ## Positive as baseline


  ## `Not taken' as baseline
  ### is this always a present level?

  data$IGRA <- as.factor(data$IGRA)
  data$PET <- as.factor(data$PET)
  data$imaging <- as.factor(data$imaging)
  data$TSTres <- as.factor(data$TSTres)
  data$IGRAorTST <- as.factor(data$IGRAorTST)

  data <- addLevel_Nottaken(data, c("IGRA","imaging","Smear","TBcult","imaging","CXR","CT","CSF","BAL","QFN","TSPOT",
                                    "HistBiop","NeedleAsp","PCR","MRI","TSTcut","PET","TSTres","IGRAorTST"))

  data$Dosanjh <- as.factor(data$Dosanjh)
  data$Dosanjh <- relevel(data$Dosanjh, ref="1")

  data$Ethnclass <- relevel(data$Ethnclass, ref="White")

  data$Country <- as.factor(data$Country)
  data$Country <- relevel(data$Country, ref="ENGLAND")

  data$Sex <- as.factor(data$Sex)


  ##TODO##
  ## what to do about TSTcut NA level?

  data[,test.names][is.na(data[,test.names]) | data[,test.names]==""] <- "Not taken"


  data$Country <- joinLevels(data$Country, list("UK"=c("ENGLAND","IRELAND","UNITED KINGDOM","WALES","SCOTLAND")))
  data$Country[data$Country=="N/A"] <- NA
  data$Country <- droplevels(data$Country)


  #  new fields -------------------------------------------------------------

  data$Alt_diag <- !data$TBconfirmed & !data$Diagoutcome%in%c("Indeterminate","Not given",NA,"")

  data$jobrisk <- data$New_occupation=="Healthcare worker"

  ## group diagnosis outcomes
  ##TODO## do over-lapping groups

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


  ## group Dosanjh category 4
  lookuplist <- list("1"=1,"2"=2,"3"=3,"4"=c("4A","4B","4C","4D"))
  data$DosanjhGrouped <- data$Dosanjh
  lfac <- levels(data$data$DosanjhGrouped)
  othrlevs <- lfac[!lfac %in% unlist(lookuplist)]
  x <- c(lookuplist, othrlevs)
  names(x) <- c(names(lookuplist), othrlevs)
  levels(data$DosanjhGrouped) <- x


  # symptoms counts ----------------------------------------------------------

  symptoms.names <- c("Cough","Fever","Ngtsweat","Wghtloss","Haemop","Leth","OtherAE1","OtherAE2","OtherAE3","OtherAE4","OtherAE5","OtherAE6")
  data$numSymptoms <- apply(data[,symptoms.names],1,sum)#, na.rm=T)


  ## group number of symptoms
  lookuplist <- list("1"=1,"2"=2,"3"=3,"4"=4, "5"=5, ">5"=c(6,7,8,9,10,11,12))
  data$numSymptomsGrouped <- as.factor(data$numSymptoms)
  lfac <- levels(data$data$numSymptomsGrouped)
  othrlevs <- lfac[!lfac %in% unlist(lookuplist)]
  x <- c(lookuplist, othrlevs)
  names(x) <- c(names(lookuplist), othrlevs)
  levels(data$numSymptomsGrouped) <- x


  symptomsEndDates.names <- paste(symptoms.names, "End", sep="")


  for (i in symptomsEndDates.names){
    data[,i] <- as.Date.POSIX(data[,i])}


  # symptoms resolved estimates ---------------------------------------------

  data$Sx_resolved <- FALSE

  ## using `all' is conservative estimates

  for (i in 1:nrow(data)){

    cols <- symptomsEndDates.names[unlist(data[i,symptoms.names])]

    if(all(is.na(cols))){   #no symptoms
      data$Sx_resolved[i] <- FALSE
    }else{
      observedEndDates <- data[i, cols]
      observedEndDates[is.na(observedEndDates)] <- as.Date("2020/01/01")       #replace missing with large dates

      if(data$step1Diag[i]){
        data$Sx_resolved[i] <- all(sapply(observedEndDates, function(x) (x - data$TBDrugStart.min[i])<drugReviewPeriod))      #response to treatment
        ## could alternatively have same as below?
        ## ...
      }else{
        # data$Sx_resolved[i] <- all(sapply(observedEndDates, function(x) x <= data[i,"TBculttestDate.orig"]))    #self-cured
        data$Sx_resolved[i] <- all(sapply(observedEndDates, function(x) x <= data[i,"TBculttestDate"]))
      }
    }
  }



  # identify out-of-sequence tests ------------------------------------------

  ## which tests are after starting treatment?
  x <- t(apply(data, 1, function(z) z[testDate.names[-1]] > z["TBDrugStart.min"]))
  colnames(x) <- tests.postDrugStart <- paste(colnames(x),".postDrugStart", sep="")

  ## which tests are after culture result?
  y <- t(apply(data, 1, function(z) z[testDate.names[-1]] > z["TBculttestDate"]))
  colnames(y) <- tests.postCult <- paste(colnames(y),".postCult", sep="")

  data <- data.frame(data, x, y)


  data <- rmPatientsInCleaning(data)

  data <- calcRiskFactorScore(data)

  data <- droplevels(data)

  return(data)
}


