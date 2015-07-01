#' Read from Access IDEA database extract.
#'
#' \code{readFromDatabaseExtract} takes Access database and updates one of the tables
#'
#' @param Patients single table to update the data total database
#' @return updated full database

readFromDatabaseExtract <- function(){

  ## read from Excel _COMPLETE_WORKBOOK_ --------------------------------------------------------

  extractNames <- readLines(system.file("extdata", "extractNames.csv", package = "IDEAdectree"))
  Patients <- read.csv(system.file("extdata", "patientdata_new_Patients_table.csv", package = "IDEAdectree"))
  patientdata_orig <- loadWorkbook(system.file("extdata", "patientdata_orig.xlsx", package = "IDEAdectree"))

  tab.names <- getSheets(patientdata_orig)
  TBlist <- sapply(1:length(tab.names), function(x) readWorksheet(patientdata_orig, sheet=x, header=TRUE, check.names=FALSE))
  names(TBlist) <- tab.names
  # y <- sapply(TBlist, function(x) x[,intersect(names(x), extractNames])     #do subset of columns before join


  ## extract only columns of interest.
  ## duplicate IDs though for multiple tb treatment entries
  # data <- join_all(TBlist[c("Patients","tb symptoms","tb treatment")], by="PatientStudyID", match="all")
  #same as:
  # data <- join_all(TBlist, by="PatientStudyID", match="all")
  # data <- unique(data[,extractNames])

  ## loss of information but no repeat IDs
  data <- join_all(TBlist, by="PatientStudyID", match="first")
  data <- data[,extractNames]

  ## update individual Excel _SHEET_ -----------------------------------------------------------

  Patients$PatientStudyID <- toupper(Patients$PatientStudyID)
  data$PatientStudyID <- toupper(data$PatientStudyID)

  patientNames <- intersect(names(Patients), extractNames)
  names <- c("PatientStudyID", setdiff(extractNames, patientNames))
  data <- merge(data[,names], Patients[,patientNames], by="PatientStudyID")

  data
}

