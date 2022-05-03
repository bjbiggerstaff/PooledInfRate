\name{summary.ipooledBin}
\alias{summary.ipooledBin}
\alias{print.summary.ipooledBin}
\alias{[.summary.ipooledBin}

\title{
Summarize ipooledBin or ipooledBinList objects
}

\description{
The functions provide more detail underlying the results from ipooledBin
}

\usage{
\method{summary}{ipooledBin}(object, ...)

\method{print}{summary.ipooledBin}(x, ...)

\method{[}{summary.ipooledBin}(x, i, j, drop = if (missing(i)) TRUE else length(cols) ==  1)

}

\arguments{
  \item{object}{an object of class \code{ipooledBin}}
  \item{x}{an object of class \code{summary.ipooledBin}}
  \item{\dots}{for additional arguments}
    \item{i,j}{elements to extract or replace. For [ and [[, these are numeric or character or, for [ only, empty or logical. Numeric values are coerced to integer as if by as.integer.}

  \item{drop}{logical. If TRUE the result is coerced to the lowest possible dimension. The default is to drop if only one column is left, but not to drop if only one row is left.}
}

\details{
Summaries include, in addition to the original object, as attributes: the number of individuals, number of pools, number of positive pools, sensitivity, and specificity, along the point estimation method used.
}

\value{
The original \code{pooledBin} object with additional list elements

\item{scale}{scale of estimates, with default 1}
\item{\dots}{future arguments}

and attributes \code{class} and \code{call}
}
\author{
Brad Biggerstaff
}

\seealso{
See also \code{\link{ipooledBin}}
}

\examples{

# Consider an imaginary example, where pools of size
# 1, 5, 10 and 50 are tested, 5 pools of each size.
# Among the 5 pools with sizes 1 and 5, no pool is positive,
# while among the 5 pools of sizes 10 and 50, 1 and 2 positive
# pools are identified, respectively.

x1 <- c(0,0,1,2)
m1 <- c(1,5,10,50)
n1 <- c(5,5,5,5)

ipooledBin(x=x1, m=m1, n=n1, sens=0.95, spec=0.99)
ipooledBin(x=x1, m=m1, n=n1, sens=0.95, spec=0.99, scale = 1000)

summary(ipooledBin(x=x1, m=m1, n=n1, sens=0.95, spec=0.99, scale=1000))

}
