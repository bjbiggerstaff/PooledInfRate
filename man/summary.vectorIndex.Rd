\name{summary.vectorIndex}
\alias{summary.vectorIndex}
\alias{summary.VI}
\alias{print.summary.vectorIndex}
\alias{print.summary.VI}
\alias{[.summary.vectorIndex}
\alias{[.summary.VI}

\title{
Summarize a vectorIndex/VI object
}
\description{
The functions provide more detail underlying the results from vectorIndex/VI
}
\usage{
\method{summary}{vectorIndex}(object, ...)
\method{summary}{VI}(object, ...)

\method{print}{summary.vectorIndex}(x, ...)
\method{print}{summary.VI}(x, ...)
%\method{print}{summary.pooledBinList}(x, ...)
%\method{print}{summary.pIRList}(x, ...)

\method{[}{summary.vectorIndex}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)
\method{[}{summary.VI}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)

}

\arguments{
  \item{object}{an object of class \code{vectorIndex}/\code{VI}}
  \item{\dots}{for additional arguments}
  \item{x}{an object of class \code{summary.vectorIndex}/\code{summary.VI}}
  \item{i,j}{elements to extract or replace. For [ and [[, these are numeric or character or, for [ only, empty or logical. Numeric values are coerced to integer as if by as.integer.}

  \item{drop}{logical. If TRUE the result is coerced to the lowest possible dimension. The default is to drop if only one column is left, but not to drop if only one row is left.}
}

\details{
In addition to the estimates, details including the vector indices for the individual species by group.

The subsetting or extractor functions \code{[.vectorIndex} and \code{[.VI} mimic the same \code{[} behavior as with data frames.

}

\value{
The original \code{vectorIndex} or \code{VI} object with additional list elements

\item{\dots}{future arguments}

and attributes \code{class} and \code{call}.
}

\author{
Brad Biggerstaff
}

\seealso{
See also \code{\link{pooledBin}}, \code{\link{pIR}}
}
\examples{
x <- c(1,0,0,0)
m <- c(50,25,10,5)

pb <- pooledBin(x,m)
summary(pb)

}
