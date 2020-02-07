
#' Clean IDEA Study Clinical Data
#'
#' A high-level function to clean IDEA clinical data for analysis.
#' There are rather a lot of assumptions and imposed contraints on values.
#' So is a bit subjective.
#'
#' @param save_data
#' @param data Individual patient records from IDEA
#'
#' @return data

cleanData <- function(data,
                      save_data = FALSE){

  # load(file=data)   #data: file address (string)

  require(lubridate)
  require(mondate)

  # TB tests ------------------------------------------------------------------

  names(data) <- gsub("testRes", "", names(data))

  ## nb culture included
  test.names.ORIGINAL <- c("QFN", "TSPOT", "TST", "Smear", "TBcult", "BAL", "HistBiop",
                           "NeedleAsp", "PCR", "CXR", "CT", "MRI")
  testDate.names.ORIGINAL <- c("QFNtestDate", "TSPOTtestDate", "TSTtestDate", "SmeartestDate",
                               "TBculttestDate", "BALtestDate", "HistBioptestDate",
                               "NeedleAsptestDate", "PCRtestDate", "CXRtestDate", "CTtestDate", "MRItestDate")


  # Discretised groups --------------------------------------------------------
  data$TSTcut <- cut(data$TST, breaks = c(0,6,15,100), right = FALSE)  #used by Yemesi @ Bham

  test.names <- c("QFN", "TSPOT", "TST", "TSTcut", "Smear", "TBcult", "CSF", "BAL",
                  "HistBiop", "NeedleAsp", "PCR", "CXR", "CT", "MRI")

  data[,c(test.names,"Othertest1Res","Othertest2Res")] <-
    data.frame(apply(data[ ,c(test.names,"Othertest1Res","Othertest2Res")], 2, toupper))

  data$IGRA <- combineTestResults(data$QFN, data$TSPOT)

  data$imaging <- combineTestResults(data$CT, data$MRI)


  # Combine with BCG status -------------------------------------------------------------

  data$TSTres <- ""
  data$TSTres[!is.na(data$TST) & data$TSTcut != "[0,6)" & data$PrevBCG == FALSE] <- "POSITIVE"
  data$TSTres[!is.na(data$TST) & data$TSTcut == "[15,100)" & data$PrevBCG == TRUE] <- "POSITIVE"
  data$TSTres[!is.na(data$TST) & data$TSTcut == "[0,6)"] <- "NEGATIVE"
  data$TSTres[!is.na(data$TST) & data$TSTcut == "[6,15)" & data$PrevBCG == TRUE] <- "NEGATIVE"

  data$IGRA[data$IGRA == "BORDERLINE"] <- NA
  data$IGRA[data$IGRA == "INDETERMINATE"] <- NA

  data$IGRAorTST <- combineTestResults(data$TSTres, data$IGRA)


  # look for PET in Othertests---------------------------------------------

  PETtest1taken <- grepl("PET", data$Othertest1Name, ignore.case = TRUE)
  PETtest2taken <- grepl("PET", data$Othertest2Name, ignore.case = TRUE)

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

  threeWeeks <- 21L
  sixWeeks   <- 42L
  twoMonths  <- 63L

  TBdrugStart.names <- c("TBdrug1Start", "TBdrug2Start", "TBdrug3Start", "TBdrug4Start", "TBdrug5Start", "TBdrug6Start")
  TBDrugEnd.names <- c("TBDrug1End", "TBDrug2End", "TBDrug3End", "TBDrug4End", "TBDrug5End", "TBDrug6End")
  testDate.names <- c("TBculttestDate","QFNtestDate", "TSPOTtestDate", "TSTtestDate", "SmeartestDate", "BALtestDate", "HistBioptestDate",
                      "NeedleAsptestDate", "PCRtestDate", "CXRtestDate", "CTtestDate", "MRItestDate", "Othertest1Date", "Othertest2Date","PETtestDate")

  for (i in c(TBdrugStart.names, TBDrugEnd.names, testDate.names,
              "EntryUK", "DateDiagCon", "DateConsent")){
    data[,i] <- as.Date.POSIX(data[,i])
  }

  ## missing dates
  data[data == "1900-01-01"] <- NA

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
      data[data$PatientStudyID == id & !is.na(data$PatientStudyID), i] <-
        min(data[data$PatientStudyID == id & !is.na(data$PatientStudyID), i], na.rm = TRUE)
    }
    for (j in TBDrugEnd.names){
      data[data$PatientStudyID == id & !is.na(data$PatientStudyID), j] <-
        min(data[data$PatientStudyID == id & !is.na(data$PatientStudyID), j], na.rm = TRUE)
    }
  }

  data <- rm.TooEarly(data, TBdrugStart.names, maxtime = threeWeeks)
  data <- rm.TooEarly(data, TBDrugEnd.names, maxtime = threeWeeks)

  data$TBDrugStart.min <- apply(data[ ,TBdrugStart.names], 1, min, na.rm = TRUE)
  data$TBDrugEnd.max   <- apply(data[ ,TBDrugEnd.names], 1, max, na.rm = TRUE)

  data <- rm.TooEarly(data, testDate.names, maxtime = threeWeeks)
  data <- calc.testDate.min(data, testDate.names, maxtime = threeWeeks)

  data$TBDrugStart.min <- as.Date.POSIX(data$TBDrugStart.min)
  data$TBDrugEnd.max   <- as.Date.POSIX(data$TBDrugEnd.max)
  data$testDate.min    <- as.Date.POSIX(data$testDate.min)

  data$DateVisitFU <- as.Date.POSIX(data$DateVisitFU)


  # lower limits culture report dates -----------------------------------------
  ## the culture takes at least _dur_ days to produce a result

  dur <- c(POSITIVE = 1,
           NEGATIVE = sixWeeks)   #days
  data$TBculttestDate.orig <- data$TBculttestDate
  dur.each <- as.vector(dur[as.character(data$TBcult)])
  data <- transform(data,
                    TBculttestDate.resMin = testDate.min + dur.each)
  data$TBculttestDate <- apply(data, 1,
                               function(x) max(x["TBculttestDate"],
                                               x["TBculttestDate.resMin"], na.rm = TRUE))
  data$TBcultCens <- data$TBculttestDate != data$TBculttestDate.orig


  # time-to-events ----------------------------------------------------------

  calcTimeToEvent <- function(testDate)
    difftime(testDate, data$testDate.min, units = "days")

  data <- transform(data,
                    TBconfirmed = (Diagoutcome %in% c("Active TB", "Active TB;Other")),

                    preCultDrug = (TBculttestDate > TBDrugStart.min),
                    preTestDrug = (testDate.min > TBDrugStart.min),

                    DrugCult_diff = round(difftime(TBculttestDate, TBDrugStart.min, units = "days")),
                    testDrug_diff = calcTimeToEvent(TBDrugStart.min),
                    testDiagCon_diff = round(calcTimeToEvent(DateDiagCon)),
                    testCult_diff = calcTimeToEvent(TBculttestDate),
                    testCultorig_diff = calcTimeToEvent(TBculttestDate.orig),
                    EntryUKtest_diff = year(testDate.min) - EntryUK_year,

                    start.to.Drug = calcTimeToEvent(TBDrugStart.min), #duplication but clearer later-on
                    start.to.Smear = calcTimeToEvent(SmeartestDate),
                    start.to.HistBiop = calcTimeToEvent(HistBioptestDate),
                    start.to.NeedleAsp = calcTimeToEvent(NeedleAsptestDate),
                    start.to.TBcultorig = calcTimeToEvent(TBculttestDate.orig),
                    start.to.TBcult = calcTimeToEvent(TBculttestDate),
                    start.to.BAL = calcTimeToEvent(BALtestDate),
                    start.to.PCR = calcTimeToEvent(PCRtestDate),
                    start.to.CXR = calcTimeToEvent(CXRtestDate),
                    start.to.TST = calcTimeToEvent(TSTtestDate),
                    start.to.CT = calcTimeToEvent(CTtestDate),
                    start.to.MRI = calcTimeToEvent(MRItestDate),
                    start.to.PET = calcTimeToEvent(PETtestDate),
                    start.to.QFN = calcTimeToEvent(QFNtestDate),
                    start.to.TSPOT = calcTimeToEvent(TSPOTtestDate),

                    start.to.other1 = calcTimeToEvent(Othertest1Date),
                    start.to.other2 = calcTimeToEvent(Othertest2Date),

                    start.to.FU = calcTimeToEvent(DateVisitFU)
  )

  # View(data.frame(data$PatientStudyID, data$testDate.min, data$DateVisitFU, data$DateVisitFU0, data$start.to.FU)[order(data$PatientStudyID),])  #check

  data$start.to.Imaging <- pmax(data$start.to.CT, data$start.to.CXR, data$start.to.MRI, data$start.to.PET, na.rm = TRUE)
  data$start.to.IGRA  <- pmax(data$start.to.QFN, data$start.to.TSPOT, na.rm = TRUE)
  data$start.to.other <- pmax(data$start.to.other1, data$start.to.other2, na.rm = TRUE)

  #NeedleAsp isnt a test but a sampling technique so could be used for other test
  data$start.to.Histology <- pmax(data$start.to.HistBiop, data$start.to.NeedleAsp, na.rm = TRUE)

  data$start.to.clinicalfeatures <- with(data, pmax(start.to.TBcultorig, start.to.Smear, start.to.Histology,
                                                    start.to.BAL, start.to.PCR, start.to.TST, start.to.Imaging,
                                                    start.to.IGRA, start.to.other, na.rm = TRUE))

  ## only interested in 2 month followup
  # data$start.to.FU[data$VisitFU!="2 month FU"] <- NA

  data$preTestDrug[is.na(data$preTestDrug)] <- FALSE

  data$step1Diag <- (data$preCultDrug & !data$preTestDrug)
  data$step1Diag[is.na(data$step1Diag)] <- FALSE


  ## remove previous TB and LTBI patients
  ##TODO##
  ## whats the field for this??

  data <- fillInEndOfTreatmentDate(data)

  data$TBDrug_diff <- difftime(data$TBDrugEnd.max, data$TBDrugStart.min, units = "days")

  ## this isn't perfect because there may be treatment gaps but is an ok approximation
  drugReviewPeriod <- twoMonths
  data$treatResponse <- (as.numeric(data$TBDrug_diff) > drugReviewPeriod)
  data$testDrug_diff_plus63days <- data$testDrug_diff + drugReviewPeriod

  TBdrugStart.freq <- getDateFrequencies(TBdrugStart.names, data)
  TBDrugEnd.freq <- getDateFrequencies(TBDrugEnd.names, data)
  testDate.freq  <- getDateFrequencies(testDate.names, data)

  ## should we fill-in missing dates with 2 month estimates?
  # data$start.to.FUmissing <- is.na(data$start.to.FU)
  # data$start.to.FU <- ifelse(data$start.to.FUmissing, data$testDrug_diff_plus63days, data$start.to.FU)

  data <- estimateTimeToDiagnosis(data)

  ##TODO##
  ## only include tests that report _before_ culture reports
  ## if they're after then substitute in "Not taken"
  # names(data)[grep("testDate", names(data))]


  ## relevel factors -----------------------------------------------------------------------
  ## Positive as baseline


  ## `Not taken' as baseline
  ### is this always a present level?

  data$IGRA <- as.factor(data$IGRA)
  data$PET  <- as.factor(data$PET)
  data$imaging <- as.factor(data$imaging)
  data$TSTres  <- as.factor(data$TSTres)
  data$IGRAorTST <- as.factor(data$IGRAorTST)

  data <- addLevel_Nottaken(data, c("IGRA","imaging","Smear","TBcult","imaging","CXR","CT","CSF","BAL","QFN","TSPOT",
                                    "HistBiop","NeedleAsp","PCR","MRI","TSTcut","PET","TSTres","IGRAorTST"))

  data$Dosanjh <- as.factor(data$Dosanjh)
  data$Dosanjh <- relevel(data$Dosanjh, ref = "1")

  data$Ethnclass <- relevel(data$Ethnclass, ref = "White")

  data$Country <- as.factor(data$Country)
  data$Country <- relevel(data$Country, ref = "ENGLAND")

  data$Sex <- as.factor(data$Sex)

  ##TODO##
  ## what to do about TSTcut NA level?

  data[ ,test.names.ORIGINAL][is.empty(data[ ,test.names.ORIGINAL]) & !is.na(data[ ,testDate.names.ORIGINAL])] <- "INDETERMINATE"
  data[ ,test.names][is.empty(data[ ,test.names])] <- "Not taken"

  data$Country <- joinLevels(data$Country,
                             list("UK" = c("ENGLAND","IRELAND","UNITED KINGDOM","WALES","SCOTLAND")))
  data$Country[data$Country == "N/A"] <- NA
  data$Country <- droplevels(data$Country)


  #  new fields -------------------------------------------------------------

  data$Alt_diag <- !data$TBconfirmed & !data$Diagoutcome %in% c("Indeterminate","Not given",NA,"")

  data$jobrisk <- data$New_occupation == "Healthcare worker"

  data <- GroupDiagOutcomes(data)

  ## group Dosanjh category 4
  lookuplist <- list("1"=1, "2"=2, "3"=3, "4"=c("4A","4B","4C","4D"))
  data$DosanjhGrouped <- data$Dosanjh
  lfac <- levels(data$data$DosanjhGrouped)
  othrlevs <- lfac[!lfac %in% unlist(lookuplist)]
  x <- c(lookuplist, othrlevs)
  names(x) <- c(names(lookuplist), othrlevs)
  levels(data$DosanjhGrouped) <- x


  # symptoms counts ----------------------------------------------------------

  symptoms.names <- c("Cough","Fever","Ngtsweat","Wghtloss","Haemop","Leth",
                      "OtherAE1","OtherAE2","OtherAE3","OtherAE4","OtherAE5","OtherAE6")
  data$numSymptoms <- apply(data[ ,symptoms.names], 1, sum)#, na.rm=T)

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
    data[,i] <- as.Date.POSIX(data[ ,i])}


  data <- EstimateSymptomsResolved(data, symptomsEndDates.names, symptoms.names, drugReviewPeriod)

  data <- is.diagByClinicalFeatures(data)


  # estimate counts of each test taken --------------------------------------

  tests <- c("TBcult", "Smear", "IGRA", "TSTres", "TSPOT", "QFN", "CXR",
             "CSF", "BAL", "HistBiop", "NeedleAsp", "PCR", "CT", "MRI", "PET")

  NumEachTest <- data.frame(apply(data[,tests] != "Not taken", 2, as.numeric))
  NumEachTest[is.na(NumEachTest)] <- 0
  names(NumEachTest) <- paste("Num_", tests, sep = "")
  data <- data.frame(data, NumEachTest)


  ## this assumes all test dates are recorded, which they aren't!
#   # identify out-of-sequence tests ------------------------------------------
#
#   ## which tests are after starting treatment?
#   x <- t(apply(data, 1, function(z) z[testDate.names[-1]] > z["TBDrugStart.min"]))
#   colnames(x) <- tests.postDrugStart <- paste(colnames(x),".postDrugStart", sep="")
#
#   ## which tests are after culture result?
#   y <- t(apply(data, 1, function(z) z[testDate.names[-1]] > z["TBculttestDate"]))
#   colnames(y) <- tests.postCult <- paste(colnames(y),".postCult", sep="")
#   data <- data.frame(data, x, y)

  data <- rmPatientsInCleaning(data)

  data <- calcRiskFactorScore(data)

  ##TODO##
  #where has this variable gone?
  # data$PatientWgt[data$PatientWgt==0] <- NA

  data <- tidyBALSputumdata(data)
  data <- tidyCXRdata(data)
  data <- tidyCTdata(data)

  ##TODO##
  #COSTS argument is now from a distribution done later in a sensitivity analysis
  # data <- data.frame(data, totalcost=calcPatientCostofTests(data))

  data <- droplevels(data)

  if (save_data){

  }

  return(data)
}


