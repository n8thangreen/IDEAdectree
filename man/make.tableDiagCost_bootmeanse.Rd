% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enhanceddiagcost_table.R
\name{make.tableDiagCost_bootmeanse}
\alias{make.tableDiagCost_bootmeanse}
\title{Make a table of time to diagnosis and total test costs}
\usage{
make.tableDiagCost_bootmeanse(data,
  sampleFUN = calc.enhancedDiagTimeAndCost.bootsample,
  iterations = 1000)
}
\arguments{
\item{sampleFUN}{Function name to subsample, or bootstrap (default) data; could use \code{identity} function}

\item{interations}{}
}
\value{
Matrix of means and sd by Dosanjh category
}
\description{
Bootstrap patient data for mean and standard deviation estimates
}
