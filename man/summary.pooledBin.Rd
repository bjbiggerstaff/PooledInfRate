\name{summary.pooledBin}
\alias{summary.pooledBin}
\alias{summary.pIR}
\alias{print.summary.pooledBin}
\alias{print.summary.pIR}
\alias{[.summary.pooledBin}
\alias{[.summary.pIR}

\title{
Summarize a pooledBin/pIR or pooledBinList/pIRList objects
}
\description{
The functions provide more detail underlying the results from pooledBin/pIR
}
\usage{
\method{summary}{pooledBin}(object, simple = FALSE, ...)
\method{summary}{pIR}(object, simple = FALSE, ...)

\method{print}{summary.pooledBin}(x, ...)
\method{print}{summary.pIR}(x, ...)
%\method{print}{summary.pooledBinList}(x, ...)
%\method{print}{summary.pIRList}(x, ...)

\method{[}{summary.pooledBin}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)
\method{[}{summary.pIR}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)

}

\arguments{
  \item{object}{an object of class \code{pooledBin}/\code{pIR}}
  \item{\dots}{for additional arguments}
  \item{x}{an object of class \code{summary.pooledBin}/\code{summary.pIR}}
  \item{simple}{logical.  If TRUE return a data frame with estimates, confidence intervals, total number of individuals, number of positive pools, and scale, by group; if FALSE return a more complete description of the call.}
  \item{i,j}{elements to extract or replace. For [ and [[, these are numeric or character or, for [ only, empty or logical. Numeric values are coerced to integer as if by as.integer.}

  \item{drop}{logical. If TRUE the result is coerced to the lowest possible dimension. The default is to drop if only one column is left, but not to drop if only one row is left.}
}

\details{
In addition to the estimates and confidence intervals, details including the total number of individuals included, the number of pools, and the number of "positive" pools are given.  Printing the result reports these, along with which point and confidence interval methods were used and the call.

The subsetting or extractor functions \code{[.pooledBin} and \code{[.pIR} mimic the same \code{[} behavior as with data frames.

}

\value{
The original \code{pooledBin} or \code{pIR} object with additional list elements

\item{scale}{scale of estimates, with default 1}
\item{pt.method}{the specified point estimate method to be used}
\item{ci.method}{the specified confidence interval methods to be used}
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
