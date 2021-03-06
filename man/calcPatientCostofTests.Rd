% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcPatientCostofTests.R
\name{calcPatientCostofTests}
\alias{calcPatientCostofTests}
\title{Sum the Cost of All Specified Tests for Each Patient}
\usage{
calcPatientCostofTests(data, COSTS, x3 = TRUE)
}
\arguments{
\item{data}{IDEA study data set}

\item{COSTS}{Named vector each test/procedure and unit costs.
If a list of cost is provided then the named mean value is used.}

\item{x3}{Should the smear and culture be multiplied by 3?
It is standard practice to perform 3 of each with a sputum sample.}
}
\value{
Single value per patient (vector)
}
\description{
\sum freq_i \times count_i \times cost_i
}
\details{
$freq$ are when a single entry in the dataset represent more that one unit cost.
$count$ are the numbeer of separate occasion at test is taken.
$cost$ is the unit cost.

Most frequencies are singular but some or >1.
Some costs consist of more than element too and these are hard-coded in to the function;
If a BAL is taken we include the cost of the procedure and the additional
cost of a culture and smear if used for these.
}
