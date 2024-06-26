% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{samplev_cpp}
\alias{samplev_cpp}
\title{Calculate a Sample Vector}
\usage{
samplev_cpp(probs, m, seed = 1L)
}
\arguments{
\item{probs}{A numeric vector of probabilities.}

\item{m}{An integer indicating the number of samples to return.}

\item{seed}{An optional integer that sets the random seed for reproducibility.}
}
\value{
A numeric vector of sampled values.
}
\description{
This function takes a vector of probabilities and an integer value \code{m}, and optionally a seed,
and returns a sampled vector based on the provided probabilities.
}
\examples{
samplev_cpp(c(0.1, 0.9), 5)
}
