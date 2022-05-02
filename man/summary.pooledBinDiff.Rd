\name{summary.pooledBinDiff}
\alias{summary.pooledBinDiff}
\alias{summary.pIRDiff}
\alias{print.summary.pooledBinDiff}
\alias{[.summary.pooledBinDiff}
\alias{print.summary.pIRDiff}
\alias{[.summary.pIRDiff}

\title{
Summarize \code{pooledBinDiff}/\code{pIRDiff} objects
}
\description{
The functions provide more detail underlying the results from pooledBinDiff/pIRDIff
}
\usage{
\method{summary}{pooledBinDiff}(object, scale = attr(x, "scale"), ...)
\method{summary}{pIRDiff}(object, scale = attr(x, "scale"), ...)

\method{print}{summary.pooledBinDiff}(x, \dots)
\method{[}{summary.pooledBinDiff}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)
\method{print}{summary.pIRDiff}(x, \dots)
\method{[}{summary.pIRDiff}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)
}


\arguments{
  \item{object}{object of class \code{pooledBinDiff} or \code{pIRDiff} }
  \item{scale}{the scale to report the results in, with default the scale passed by the object}
  \item{x}{object of class \code{summary.pooledBinDiff}/\code{summary.pIRDiff}}
    \item{i,j}{elements to extract or replace. For [ and [[, these are numeric or character or, for [ only, empty or logical. Numeric values are coerced to integer as if by as.integer.}

  \item{drop}{logical. If TRUE the result is coerced to the lowest possible dimension. The default is to drop if only one column is left, but not to drop if only one row is left.}
  \item{\dots}{future arguments}
}


\details{
In addition to the estimates and confidence intervals for the difference, details for the indivdiual gropus including the total number of individuals included, the number of pools, and the number of "positive" pools are given.  Printing the result reports these, along with which point and confidence interval methods were used and the call.

The subsetting or extractor functions \code{[.pooledBin} and \code{[.pIR} mimic the same \code{[} behavior as with data frames.

}

\value{
The original \code{pooledBinDiff}/\code{pIRDiff} object with additional list elements

\item{pt.method}{the specified point estimate method to be used}
\item{ci.method}{the specified confidence interval methods to be used}

and attributes \code{class}, \code{grp.pooledBin}, \code{PtEstName}, \code{CIEstName}, and \code{scale}
}
\author{
Brad Biggerstaff
}
\seealso{
\code{\link{pooledBinDiff}}, \code{\link{pIRDiff}}
}
\examples{

x1 <- c(0,0,1,2)
m1 <- c(1,5,10,50)
n1 <- c(5,5,5,5)

x2<- c(0,1,0,4)
m2 <- m1
n2 <- n1

ex.dat <- data.frame(NumPos = c(x1,x2),
                      NumPos2 = c(1,0,0,1,1,0,0,0),
                      PoolSize = c(m1,m2),
                      NumPools = c(n1,n2),
                      Population = rep(1:2, each=4))

pbd <- pooledBinDiff(NumPos ~ PoolSize + n(NumPools) | Population, data = ex.dat)
summary(pbd)


}
