#' Clean IDEA clinical data.
#'
#' \code{cleanData} cleans IDEA clinical data for analysis
#'
#' @param data Individual patient records
#' @return data

cleanData <- function(data){
  ## data: file address (string)
  ##

  load(file=data)


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

  for (i in c(TBdrugStart.names, TBDrugEnd.names, testDate.names, "EntryUK","DateDiagCon","DateConsent")){
    data[,i] <- as.Date.POSIX(data[,i])}

  ## missing dates
  data[data=="1900-01-01"] <- NA

  data$EntryUK_year <- year(data$EntryUK)

  ## how can we know which are filled-in dates and which are the accurate ones?
  ## take the hit of omitting info for true dates too
  ## not so bad for day-month but for day this is not good
  ### ie all Junes or all 1st of month dates
  for (i in 1:ncol(data)){
    data[grepl("-06-15", data[,i]),i] <- NA}

  ## combine drug start dates for multiple patient records
  for (id in unique(data$PatientStudyID)[!is.na(unique(data$PatientStudyID))]){
    for (i in TBdrugStart.names){
      data[data$PatientStudyID==id & !is.na(data$PatientStudyID), i] <-
        min(data[data$PatientStudyID==id & !is.na(data$PatientStudyID), i], na.rm=T)
    }
  }

  for (id in unique(data$PatientStudyID)[!is.na(unique(data$PatientStudyID))]){
    for (i in TBDrugEnd.names){
      data[data$PatientStudyID==id & !is.na(data$PatientStudyID), i] <-
        min(data[data$PatientStudyID==id & !is.na(data$PatientStudyID), i], na.rm=T)
    }
  }


  data$TBDrugStart.min <- apply(data[,TBdrugStart.names], 1, min, na.rm=T)
  data$TBDrugEnd.max   <- apply(data[,TBDrugEnd.names], 1, max, na.rm=T)
  data$testDate.min    <- apply(data[,testDate.names[-1]], 1, min, na.rm=T)   #exclude culture test date in this

  data$TBDrugStart.min <- as.Date.POSIX(data$TBDrugStart.min)
  data$TBDrugEnd.max   <- as.Date.POSIX(data$TBDrugEnd.max)
  data$testDate.min    <- as.Date.POSIX(data$testDate.min)


  # lower limits culture report dates -----------------------------------------

  dur <- c(POSITIVE=1, NEGATIVE=42)   #days
  data$TBculttestDate.orig <- data$TBculttestDate
  dur.each <- as.vector(dur[as.character(data$TBcult)])
  data <- transform(data, TBculttestDate.resMin = testDate.min+dur.each)
  data$TBculttestDate <- apply(data, 1,
                               function(x) max(x["TBculttestDate"], x["TBculttestDate.resMin"], na.rm=T))
  data$TBcultCens <- data$TBculttestDate!=data$TBculttestDate.orig


  ## time-to-events
  data <- transform(data,
                    TBconfirmed = (Diagoutcome%in%c("Active TB", "Active TB;Other")),

                    preCultDrug = (TBculttestDate>TBDrugStart.min),
                    preTestDrug = (testDate.min>TBDrugStart.min),

                    DrugCult_diff = round(difftime(TBculttestDate, TBDrugStart.min, units="days")),
                    testDrug_diff = difftime(TBDrugStart.min, testDate.min, units="days"),
                    testDiagCon_diff = round(difftime(DateDiagCon, testDate.min, units="days")),
                    testCult_diff = difftime(TBculttestDate, testDate.min, units="days"),

                    EntryUKtest_diff = year(testDate.min)- EntryUK_year
  )

  data$preTestDrug[is.na(data$preTestDrug)] <- FALSE

  data$step1Diag <- (data$preCultDrug & !data$preTestDrug)
  data$step1Diag[is.na(data$step1Diag)] <- FALSE


  ## remove previous TB and LTBI patients
  ##TODO##
  ## whats the field for this??


  TBtreatdur.num <- c("6 months"=6, "2 months"=2, "3 months"=3, "9 months"=9, "12 months"=12)


  ## add-on treatment duration time when end of treatment date not given
  tdur <- is.na(data$TBDrugEnd.max) & data$TBtreatdur%in%names(TBtreatdur.num)
  data$TBDrugEnd.max[tdur] <- as.Date(as.mondate(data$TBDrugStart.min)[tdur] + TBtreatdur.num[data$TBtreatdur][tdur])
  data$TBDrugEnd.max.estimate <- tdur

  data$TBDrug_diff <- difftime(data$TBDrugEnd.max, data$TBDrugStart.min, units="days")

  ## this isn't perfect because there may be treatment gaps but is an ok approximation
  data$treatResponse <- (as.numeric(data$TBDrug_diff)>63)     #i.e. ~2 months


  TBdrugStart.freq <- getDateFrequencies(TBdrugStart.names, data)
  TBDrugEnd.freq <- getDateFrequencies(TBDrugEnd.names, data)
  testDate.freq <- getDateFrequencies(testDate.names, data)


  data$timeToDiag <- estimateTimeToDiagnosis(data)


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

  ### add level
  levels(data$IGRA) <- c(levels(data$IGRA),"Not taken")
  levels(data$imaging) <- c(levels(data$imaging),"Not taken")
  levels(data$Smear) <- c(levels(data$Smear),"Not taken")
  levels(data$TBcult) <- c( levels(data$TBcult),"Not taken")
  levels(data$imaging) <- c(levels(data$imaging),"Not taken")
  levels(data$CXR) <- c(levels(data$CXR),"Not taken")
  levels(data$CT) <- c(levels(data$CT),"Not taken")
  levels(data$CSF) <- c(levels(data$CSF),"Not taken")
  levels(data$BAL) <- c(levels(data$BAL),"Not taken")
  levels(data$QFN) <- c( levels(data$QFN),"Not taken")
  levels(data$TSPOT) <- c( levels(data$TSPOT),"Not taken")
  levels(data$HistBiop) <- c( levels(data$HistBiop),"Not taken")
  levels(data$NeedleAsp) <- c(levels(data$NeedleAsp),"Not taken")
  levels(data$PCR) <- c(levels(data$PCR),"Not taken")
  levels(data$MRI) <- c(levels(data$MRI),"Not taken")
  levels(data$TSTcut) <- c(levels(data$TSTcut),"Not taken")
  levels(data$PET) <- c(levels(data$PET),"Not taken")
  levels(data$TSTres) <- c(levels(data$TSTres),"Not taken")
  levels(data$IGRAorTST) <- c(levels(data$IGRAorTST),"Not taken")

  data$Smear <- relevel(data$Smear, ref="Not taken")
  data$TBcult <- relevel(data$TBcult, ref="Not taken")
  data$CXR <- relevel(data$CXR, ref="Not taken")
  data$CT <- relevel(data$CT, ref="Not taken")
  data$IGRA <- relevel(data$IGRA, ref="Not taken")
  data$imaging <- relevel(data$imaging, ref="Not taken")
  data$CSF <- relevel(data$CSF, ref="Not taken")
  data$BAL <- relevel(data$BAL, ref="Not taken")
  data$QFN <- relevel(data$QFN, ref="Not taken")
  data$TSPOT <- relevel(data$TSPOT, ref="Not taken")
  data$HistBiop <- relevel(data$HistBiop, ref="Not taken")
  data$NeedleAsp <- relevel(data$NeedleAsp, ref="Not taken")
  data$PCR <-relevel(data$PCR, ref="Not taken")
  data$MRI <- relevel(data$MRI, ref="Not taken")
  data$TSTcut <- relevel(data$TSTcut, ref="Not taken")
  data$PET <- relevel(data$PET, ref="Not taken")
  data$TSTres <- relevel(data$TSTres, ref="Not taken")
  data$IGRAorTST <- relevel(data$IGRAorTST, ref="Not taken")

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
  ##TODO##confirm this grouping with clinicians...
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
        data$Sx_resolved[i] <- all(sapply(observedEndDates, function(x) (x - data$TBDrugStart.min[i])<63))      #response to treatment
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


  # removals -----------------------------------------------------------------

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


  # risk factor score -------------------------------------------------------

  riskfacs <- c("WHOcut","jobrisk","numSymptoms","TBcont","PatientAge","Sex","Ethnclass","CurrHomeless","HIVpos")
  formula <- paste("TBconfirmed~", paste(riskfacs,collapse="+"), sep="")
  fit <- glm(formula, data=data, family = binomial, na.action = na.exclude)
  data$riskfacScore <- predict(fit, type = "response", na.action = na.exclude)


  data
}


