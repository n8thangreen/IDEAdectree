
#' TRUE element if NA or ""
#'
#' @param df dataframe
#'
#' @return logical array

is.empty <- function(df){
  return(is.na(df) | df=="")
}
