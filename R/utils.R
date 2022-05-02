#"bracket.root" <-
#  function(f, lower=0, upper=1, ..., scaling = 1.6, direction = c("both", "up",
#                                                                  "down"), max.iter = 50)
#  {
#    direction <- match.arg(direction)
#    if(lower == upper)
#      stop("lower cannot equal upper")
#    f1 <- f(lower, ...)
#    f2 <- f(upper, ...)
#    done <- FALSE
#    num.iter <- 0
#    while(!done) {
#      #cprint(lower)
#      #cprint(upper)
#      num.iter <- num.iter + 1
#      if(num.iter >= max.iter)
#        stop(paste("root not bracketed in", max.iter,
#                   "iterations\n"))
#      if(f1 * f2 < 0) {
#        done <- TRUE
#        ans <- c(lower, upper)
#      }
#      switch(direction,
#             both = if(abs(f1) < abs(f2)) f1 <- f(lower <- lower +
#                                                    scaling * (lower - upper), ...) else f2 <-
#               f(upper <- upper + scaling * (upper -
#                                               lower), ...),
#             down = f1 <- f(lower <- lower + scaling * (lower -
#                                                          upper), ...),
#             up = f2 <- f(upper <- upper + scaling * (upper - lower),
#                          ...))
#    }
#    sort(ans)
#  }



# Use this function to bracket a root by lower and upper
# that is bounded below by lbnd and above by ubnd
# idea is to start with (lower,upper) as brackets, then
# step toward the appropriate bound, lbnd or ubnd,
# by moving the the midpoint between lower and lbnd
# (upper and ubnd) updating lower (upper) each time
"bracket.bounded.root" <-
  function(f,lower=0,upper=1,lbnd=0,ubnd=1,...,slicing=2,max.iter=50)
  {
    if(lower == upper) stop("lower cannot equal upper")
    f1 <- f(lower,...)
    f2 <- f(upper,...)
    done <- FALSE
    num.iter <- 0
    while(!done){
      num.iter <- num.iter + 1
      if(num.iter >= max.iter) stop(paste("root not bracketed in",max.iter,"iterations\n"))
      if(f1*f2 < 0) {
        done <- TRUE
        ans <- c(lower,upper)
      }
      # march toward the bound in slicing steps
      if(abs(f1) < abs(f2))
        f1 <- f(lower <- (lower+lbnd)/slicing,...)
      else
        f2 <- f(upper <- (upper+ubnd)/slicing,...)
    }
    sort(ans)
  }


#"printc" <-
#  function(x,digits=options()$digits)
#  {
#    namex <- deparse(substitute(x))
#    cat(paste(namex,"=",round(x,digits=digits),"\n"))
#  }

# make a matrix with each row a different 'comparison' vector, e.g., c(1,0,0,-1),
# for all pairwise comparisons of n.groups groups
"allPairsMat" <- function(n.groups){
  n.pairs <- n.groups*(n.groups-1)/2
  m <- matrix(0,n.pairs,n.groups)
  xlow <- 1
  for(i in 1:(n.groups-1)){
    xind <- seq(from=xlow,length.out=(n.groups-i))
    xlow <- xind[length(xind)]+1
    yind <- seq(from=i,length.out=n.groups-i+1)
    m[xind,yind] <- cbind(1,diag(-1,n.groups-i))
  }
  m
}
