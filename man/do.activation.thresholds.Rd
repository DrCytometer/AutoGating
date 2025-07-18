% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/do_activation_thresholds.R
\name{do.activation.thresholds}
\alias{do.activation.thresholds}
\title{do.activation.thresholds}
\usage{
do.activation.thresholds(
  fcs.expr.data,
  fcs.population.gates,
  fcs.figure.dir,
  fcs.activation.thresholds,
  fcs.activation.thresholds.param
)
}
\arguments{
\item{fcs.expr.data}{The marker expression data from the fcs files.}

\item{fcs.population.gates}{List of population-defining (cell type) gates.}

\item{fcs.activation.thresholds}{The thresholds used to define the positive/negative
cutoffs for the activation marker gates.}

\item{fcs.activation.thresholds.param}{Parameters for calculating the
activation gate thresholds.}
}
\value{
description
}
\description{
A function to calculate the thresholds for gating activation markers in the panel.
}
