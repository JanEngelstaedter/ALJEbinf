% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coordinates.R
\name{getCoordinates}
\alias{getCoordinates}
\title{Calculate a table of coordinates}
\usage{
getCoordinates(dnaFocal, dnaRef, aligOutput = FALSE)
}
\arguments{
\item{dnaFocal}{Focal DNA sequence.}

\item{dnaRef}{Reference DNA sequence, e.g. Escherichia coli.}

\item{aligOutput}{If FALSE (default), only the actual coordinates are returned.
If TRUE then the alignment is also returned.}
}
\value{
By default, a table with two columns: the position in the focal sequence (posFocal),
and the position in the reference sequence (posRef).
If aligOutput==TRUE, a list containing the coordinates and the alignment.
}
\description{
This function takes two DNA sequences, aligns the translated amino acid sequences,
and calculates a table of corresponding positions in this alignment ("coordinates").
}
