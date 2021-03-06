% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyBALdata.R
\name{tidyBALSputumdata}
\alias{tidyBALSputumdata}
\title{Separate and combine BAL and sputum information}
\usage{
tidyBALSputumdata(data)
}
\arguments{
\item{data}{}
}
\value{
data
}
\description{
The BAL and sputum data are poorly recorded/structured and so
create separate smear and culture BAL fields in to 4 fields.
}
\details{
Searches on keywords e.g. "BAL" or "BRONCHIAL WASH" in \code{TypeSmear, Bodysite_smear, Othertest1Name, Othertest2Name},
and \code{TBcultsite, Bodysite_culture, Othertest1Name, Othertest2Name}.

We assume that if there was a BAL then there was also a previous negative sputum sample,
(ignore the possibility if they couldn't take a sputum sample)
and if there was a BAL there this was for both smear and culture.
}
