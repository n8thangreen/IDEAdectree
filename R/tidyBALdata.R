
#' Separate and combine BAL and sputum information
#'
#' The BAL and sputum data are poorly recorded/structured and so
#' create separate smear and culture BAL fields in to 4 fields.
#'
#' Searches on keywords e.g. "BAL" or "BRONCHIAL WASH" in \code{TypeSmear, Bodysite_smear, Othertest1Name, Othertest2Name},
#' and \code{TBcultsite, Bodysite_culture, Othertest1Name, Othertest2Name}.
#'
#' We assume that if there was a BAL then there was also a previous negative sputum sample,
#' (ignore the possibility if they couldn't take a sputum sample)
#' and if there was a BAL there this was for both smear and culture.
#'
#' @param data
#'
#' @return data

tidyBALSputumdata <- function(data){

#   sputumWords <- c("SPUTUM", "SPITUM")
#   BALWords <- c("BAL", "BRONCHIAL WASH", "BRONCO LAVAGE", "BRONCHOSCOPY")
  # SmearWords <- MICROSCOPY
#
#   data$SmearBAL <-grepl(paste(BALWords,collapse="|"), data$TypeSmear, ignore.case=T)|
#     grepl(paste(BALWords,collapse="|"), data$Bodysite_smear, ignore.case=T)

  ## which Smears are BAL?
  data$SmearBAL <- (grepl("BAL", data$TypeSmear, ignore.case=T)|
                    grepl("BRONCHIAL WASH", data$TypeSmear, ignore.case=T)|
                    grepl("BAL", data$Bodysite_smear, ignore.case=T)|
                    grepl("BRONCHIAL WASH", data$Bodysite_smear, ignore.case=T))

  data$SmearSputum <- (grepl("SPUTUM", data$TypeSmear, ignore.case=T)|
                       grepl("SPITUM", data$TypeSmear, ignore.case=T)|
                       grepl("SPUTUM", data$Bodysite_smear, ignore.case=T)|
                       grepl("SPITUM", data$Bodysite_smear, ignore.case=T))

  ## which other tests are smear BAL?
  data$SmearBAL2 <- (grepl("BAL", data$Othertest1Name, ignore.case=T)|
                     grepl("BRONCHIAL WASH", data$Othertest1Name, ignore.case=T)|
                     grepl("BRONCO LAVAGE", data$Othertest1Name, ignore.case=T)|
                     grepl("BRONCHOSCOPY", data$Othertest2Name, ignore.case=T)) &
    grepl("SMEAR", data$Othertest1Name, ignore.case=T)

  data$SmearSputum2 <- (grepl("SPUTUM", data$Othertest1Name, ignore.case=T)|
                        grepl("SPITUM", data$Othertest1Name, ignore.case=T)) &
    grepl("SMEAR", data$Othertest1Name, ignore.case=T)

  data$SmearBAL3 <- (grepl("BAL", data$Othertest2Name, ignore.case=T)|
                     grepl("BRONCHIAL WASH", data$Othertest2Name, ignore.case=T)|
                     grepl("BRONCO LAVAGE", data$Othertest2Name, ignore.case=T)|
                     grepl("BRONCHOSCOPY", data$Othertest2Name, ignore.case=T)) &
    grepl("SMEAR", data$Othertest2Name, ignore.case=T)

  data$SmearSputum3 <- (grepl("SPUTUM", data$Othertest2Name, ignore.case=T)|
                        grepl("SPITUM", data$Othertest2Name, ignore.case=T)) &
    grepl("SMEAR", data$Othertest2Name, ignore.case=T)


  ## which culture are BAL?
  data$TBcultBAL <- (grepl("BAL", data$TBcultsite, ignore.case=T)|
                       grepl("BRONCHIAL WASH", data$TBcultsite, ignore.case=T)|
                       grepl("BAL", data$Bodysite_culture, ignore.case=T)|
                       grepl("BRONCHIAL WASH", data$Bodysite_culture, ignore.case=T))

  data$TBcultSputum <- (grepl("BAL", data$TBcultsite, ignore.case=T)|
                          grepl("BRONCHIAL WASH", data$TBcultsite, ignore.case=T)|
                          grepl("BAL", data$Bodysite_culture, ignore.case=T)|
                          grepl("BRONCHIAL WASH", data$Bodysite_culture, ignore.case=T))

  ## which other tests are culture BAL?
  data$TBcultBAL2 <- (grepl("BAL", data$Othertest1Name, ignore.case=T)|
                        grepl("BRONCHIAL WASH", data$Othertest1Name, ignore.case=T)|
                        grepl("BRONCO LAVAGE", data$Othertest1Name, ignore.case=T)|
                        grepl("BRONCHOSCOPY", data$Othertest2Name, ignore.case=T)) &
    grepl("CULTURE", data$Othertest1Name, ignore.case=T)

  data$TBcultSputum2 <- (grepl("SPUTUM", data$Othertest1Name, ignore.case=T)|
                           grepl("SPITUM", data$Othertest1Name, ignore.case=T)) &
    grepl("CULTURE", data$Othertest1Name, ignore.case=T)

  data$TBcultBAL3 <- (grepl("BAL", data$Othertest2Name, ignore.case=T)|
                        grepl("BRONCHIAL WASH", data$Othertest2Name, ignore.case=T)|
                        grepl("BRONCO LAVAGE", data$Othertest2Name, ignore.case=T)|
                        grepl("BRONCHOSCOPY", data$Othertest2Name, ignore.case=T)) &
    grepl("CULTURE", data$Othertest2Name, ignore.case=T)

  data$TBcultSputum3 <- (grepl("SPUTUM", data$Othertest2Name, ignore.case=T)|
                           grepl("SPITUM", data$Othertest2Name, ignore.case=T)) &
    grepl("CULTURE", data$Othertest2Name, ignore.case=T)


  ## associated Dates
  data$SmearBALdate <- data$SmearSputumdate <- data$SmeartestDate
  data$SmearBALdate[!data$SmearBAL] <- NA
  data$SmearSputumdate[!data$SmearSputum] <- NA

  data$SmearBALdate2 <- data$SmearSputumdate2 <- data$Othertest1Date
  data$SmearBALdate2[!data$SmearBAL2] <- NA
  data$SmearSputumdate2[!data$SmearSputum2] <- NA

  data$SmearBALdate3 <- data$SmearSputumdate3 <- data$Othertest2Date
  data$SmearBALdate3[!data$SmearBAL3] <- NA
  data$SmearSputumdate3[!data$SmearSputum3] <- NA


  data$TBcultBALdate <- data$TBcultSputumdate <- data$TBculttestDate
  data$TBcultBALdate[!data$TBcultBAL] <- NA
  data$TBcultSputumdate[!data$TBcultSputum] <- NA

  data$TBcultBALdate2 <- data$TBcultSputumdate2 <- data$Othertest1Date
  data$TBcultBALdate2[!data$TBcultBAL2] <- NA
  data$TBcultSputumdate2[!data$TBcultSputum2] <- NA

  data$TBcultBALdate3 <- data$TBcultSputumdate3 <- data$Othertest2Date
  data$TBcultBALdate3[!data$TBcultBAL3] <- NA
  data$TBcultSputumdate3[!data$TBcultSputum3] <- NA


  # fill-in tests taken and results according to assumptions ---------------------

  atleastone.SmearSputum <- with(data, SmearSputum|SmearSputum2|SmearSputum3)
  atleastone.TBcultSputum <- with(data, TBcultSputum|TBcultSputum2|TBcultSputum3)
  atleastone.BAL <- with(data, SmearBAL|SmearBAL2|SmearBAL3|TBcultBAL|TBcultBAL2|TBcultBAL3)
  NoTBcultBAL <- with(data, !TBcultBAL&!TBcultBAL2&!TBcultBAL3)
  NoSmearBAL  <- with(data, !SmearBAL&!SmearBAL2&!SmearBAL3)

  data$SmearSputum[!atleastone.SmearSputum & atleastone.BAL]  <- TRUE       #no recorded smears
  data$TBcultSputum[!atleastone.TBcultSputum & atleastone.BAL] <- TRUE

  data$SmearBAL[NoSmearBAL & !NoTBcultBAL]  <- TRUE               #at least one culture BAL
  data$TBcultBAL[!NoSmearBAL & NoTBcultBAL] <- TRUE


  data$SmearSputumRes <- data$Smear
  data$SmearSputumRes[!data$SmearSputum] <- NA
  data$SmearSputumRes[!data$SmearSputum & atleastone.BAL] <- "NEGATIVE"
  data$TBcultSputumRes <- data$TBcult
  data$TBcultSputumRes[!data$TBcultSputum] <- NA
  data$TBcultSputumRes[!data$TBcultSputum & atleastone.BAL] <- "NEGATIVE"

  data$SmearBALRes <- data$Smear
  data$SmearBALRes[!data$SmearBAL] <- NA
  data$TBcultBALRes <- data$TBcult
  data$TBcultBALRes[!data$TBcultBAL] <- NA


  data$Num_TBcultSmearBAL <- atleastone.BAL


  ##TODO##
  ## what do we do with the BAL columns, because we don't know if its culture or smear?

  return(data)
}

# what is proportion of bal/sputum only?
# table(data$TBcultBAL, data$SmearBAL, toupper(data$BAL)!="NOT TAKEN")
# View(data[(data$TBcultBAL==TRUE | data$SmearBAL==TRUE) & toupper(data$BAL)=="NOT TAKEN",])


#' Tidy multiple chest X-ray data
#'
#' @return data

tidyCXRdata <- function(data){

  dict <- c("CXR", "CHEST X-RAY", "CHEST XRAY")

  data$CXR2 <- grepl(paste(dict,collapse="|"), data$Othertest1Name, ignore.case=T)
  data$CXR3 <- grepl(paste(dict,collapse="|"), data$Othertest2Name, ignore.case=T)

  data$CXR2date <- data$CXR3date <- NA
  data$CXR2date[data$CXR2] <- data$Othertest1Date[data$CXR2]
  data$CXR3date[data$CXR3] <- data$Othertest2Date[data$CXR3]

  data$Num_CXR <- data$Num_CXR + as.numeric(data$CXR2) + as.numeric(data$CXR3)
  data
}


#' Tidy multiple chest X-ray data
#'
#' @return data

tidyCTdata <- function(data){

  dict <- c("CT ", " CT")

  data$CT2 <- grepl(paste(dict,collapse="|"), data$Othertest1Name, ignore.case=T)
  data$CT3 <- grepl(paste(dict,collapse="|"), data$Othertest2Name, ignore.case=T)

  data$CT2date <- data$CT3date <- NA
  data$CT2date[data$CXR2] <- data$Othertest1Date[data$CT2]
  data$CT3date[data$CXR3] <- data$Othertest2Date[data$CT3]

  data$Num_CT <- data$Num_CT + as.numeric(data$CT2) + as.numeric(data$CT3)
  data
}


#' Combine CSF information
#'
#' The CSF data is poorly recorded/structured and so
#' create separate smear and culture CSF fields in to 2 types BAL/non-BAL.
#'
#' Searches on keywords "CSF" in TypeSmear, Bodysite_smear, Othertest1Name and Othertest2Name,
#' and TBcultsite, Bodysite_culture, Othertest1Name and Othertest2Name.
#' Used for extra-pulmonary TB only.
#'
#' TODO
#'
#' @param data
#'
#' @return data

tidyCSFdata <- function(data){

  ##TODO##
}
