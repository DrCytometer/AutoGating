% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_experiment.R
\name{read.experiment}
\alias{read.experiment}
\title{read.experiment}
\usage{
read.experiment(path, fcs.experiment.filename)
}
\arguments{
\item{path}{The file path to the location of the fcs.experiment.filename file.}

\item{fcs.experiment.filename}{The csv file defining the experiment.}
}
\value{
fcs.experiment
\describe{
 \item{experiment}{Definition of the experiment.}
 \item{sample}{Character vector of the fcs samples.}
 \item{sample.n}{Numeric. Number of samples.}
 \item{sample.filename}{Character vector of the fcs file names.}
 \item{sample.gate.calculation}{Samples to be used in calculating the gates.}
 \item{sample.gate.calc.filename}{Filenames of samples to be use for gate calculation.}
 \item{sample.gate.calculation.n}{Numeric. Number of samples for gate calculation.}
}
}
\description{
The function reads in the user's definition of their experiment from the csv
file and outputs a list with the information structured for downstream use.
}
