\name{summary.pooledBin}
\alias{summary.pooledBin}
\alias{summary.pIR}
\alias{summary.pooledBinList}
\alias{summary.pIRList}
\alias{print.summary.pooledBin}
\alias{print.summary.pooledBinList}
\alias{print.summary.pIR}
\alias{print.summary.pIRList}

\title{
Summarize a pooledBin/pIR or pooledBinList/pIRList objects
}
\description{
The functions provide more detail underlying the results from pooledBin/pIR
}
\usage{
\method{summary}{pooledBin}(object, ...)
\method{summary}{pIR}(object, ...)

\method{print}{summary.pooledBin}(x, ...)
\method{print}{summary.pIR}(x, ...)
\method{print}{summary.pooledBinList}(x, ...)
\method{print}{summary.pIRList}(x, ...)

}

\arguments{
  \item{object}{an object of class \code{pooledBin}/\code{pIR} or \code{pooledBinList}/\code{pIRList}}
  \item{\dots}{for additional arguments}
  \item{x}{an object of class \code{summary.pooledBin}/\code{summary.pIR}} or \code{summary.pooledBinList}/\code{summary.pIRList}
}

\details{
In addition to the estimates and confidence intervals, details including the total number of individuals included, the number of pools, and the number of "positive" pools are given.  Printing the result reports these, along with which point and confidence interval methods were used and the call.
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
