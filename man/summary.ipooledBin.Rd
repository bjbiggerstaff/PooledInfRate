\name{summary.ipooledBin}
\alias{summary.ipooledBin}
\alias{summary.ipooledBinList}
\alias{print.summary.ipooledBin}
\alias{print.summary.ipooledBinList}

\title{
Summarize ipooledBin or ipooledBinList objects
}

\description{
The functions provide more detail underlying the results from ipooledBin
}

\usage{
\method{summary}{ipooledBin}(object, ...)

\method{print}{summary.ipooledBin}(x, ...)
\method{print}{summary.ipooledBinList}(x, ...)

}

\arguments{
  \item{object}{an object of class \code{ipooledBin} or \code{ipooledBinList}}
  \item{\dots}{for additional arguments}
  \item{x}{an object of class \code{summary.ipooledBin} or \code{summary.ipooledBinList}}
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
