
#' Sum the Cost of All Specified Tests for Each Patient
#'
#' \deqn{\sum freq_i \times count_i \times cost_i}
#'
#' $freq$ are when a single entry in the dataset represent more that one unit cost.
#' $count$ are the numbeer of separate occasion at test is taken.
#' $cost$ is the unit cost.
#'
#' Most frequencies are singular but some or >1.
#' Some costs consist of more than element too and these are hard-coded in to the function;
#' If a BAL is taken we include the cost of the procedure and the additional
#' cost of a culture and smear if used for these.
#'
#' @param data IDEA study data set
#' @param COSTS Named vector each test/procedure and unit costs.
#'              If a list of cost is provided then the named mean value is used.
#' @param x3 Should the smear and culture be multiplied by 3?
#'           It is standard practice to perform 3 of each with a sputum sample.
#'
#' @return Single value per patient (vector)
#' @export
#'
calcPatientCostofTests <- function(data,
                                   COSTS,
                                   x3 = TRUE){

  x3 <-
    assertive.base::use_first(x3) %>%
    assertive.base::coerce_to("logical")

  freq <- ifelse(x3, 3, 1)
  if (is(COSTS, "list")) { COSTS <- COSTS[["mean"]] }

  namesNumEachTest <- grep(names(data),
                           pattern = "Num_",
                           value = TRUE)

  testcosts <- c(TBcult = freq*COSTS$TBcult,
                 Smear = freq*COSTS$Smear,
                 IGRA  = (COSTS$TSPOT + COSTS$QFN)/2,  #COSTS$IGRA,
                 TSTres = COSTS$TSTres,
                 TSPOT  = COSTS$TSPOT,
                 QFN = COSTS$QFN,
                 CXR = COSTS$CXR,
                 CSF = COSTS$CSF,
                 BAL = COSTS$BAL,
                 HistBiop  = COSTS$HistBiop,
                 NeedleAsp = COSTS$NeedleAsp,
                 PCR = COSTS$PCR,
                 CT  = COSTS$CT,
                 MRI = COSTS$MRI,
                 PET = COSTS$PET,
                 EBUS = COSTS$EBUS,
                 TBcultSmearBAL = COSTS$TBcult + COSTS$Smear)#,
  # Risk_factors=0)  #why did I include this??

  ## if no test count data
  if (length(namesNumEachTest) == 0) {

    data <- as.matrix(data) # coherse because factor

    data[data %in% c("POSITIVE", "NEGATIVE")] <- 1L  # taken test
    data[data != 1 | is.na(data)] <- 0L              # not take

    class(data) <- "numeric"

    res <- testcosts[colnames(data)] %*% t(data)

  }else{

    ## cross-product to account for multiple
    ## same tests per patient
    Num_names <- paste("Num_", names(testcosts), sep = "")
    res <- c(testcosts %*% t(data[ ,Num_names]))
  }

  return(res)
}

