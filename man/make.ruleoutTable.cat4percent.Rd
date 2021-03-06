% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_ruleoutTable-bysubgroups.R
\name{make.ruleoutTable.cat4percent}
\alias{make.ruleoutTable.cat4percent}
\title{Make rule-out table for multiple cohort prevalences}
\usage{
make.ruleoutTable.cat4percent(cat4percent.seq = seq(1, 100, by = 1),
  Ctest)
}
\arguments{
\item{Ctest}{}
}
\value{
out
}
\description{
Iterate over a range of prevalences 0-100% and row-bind results (potentially a high number of rows).
}
