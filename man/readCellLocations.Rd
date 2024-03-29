% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readFunctions.R
\name{readCellLocations}
\alias{readCellLocations}
\title{Load cell locations file}
\usage{
readCellLocations(filepath, num_dims = 3)
}
\arguments{
\item{filepath}{the path to the containing results folder}

\item{num_dims}{the number of dimensions of the simulation}
}
\value{
A list of data frames, one for each time step, containing the cell locations.
}
\description{
Read in the results.viznodes file and split into locations
}
