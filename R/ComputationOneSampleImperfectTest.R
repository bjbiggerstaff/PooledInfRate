# These functions will later be incorporated into the main functionality of pooledBin/pIR (and vectorIndex/VI)
# but because there are not corresponding CIs, I'm holding off on doing that, and just including them
# in the package namespace for direct use
###################################################################################################
###################################################################################################
#### imperfect test estimation functions ----
###################################################################################################
###################################################################################################
# Newton-Raphson iteration - so don't need root-finder
"ipIR.mle" <- "ipooledbinom.mle" <- function(x,m,n=rep(1,length(m)),sens=1,spec=1, tol=1e-12,max.iter=10000){
  # a, b match the manuscript; but sens, spec are more intuitive as function arguments
  a <- sens
  b <- spec

  # if really just the regular binomial, return the simple proportion
  if(all(m==1) & all(a == 1) & all(b == 1)) return(sum(x)/sum(n))
  if(all(a==1) & all(b==1)) return(pooledbinom.mle(x,m,n))
  if(sum(x) == 0) return(0)
  # will add a check that a solution exists here...
  r <- a + b - 1
  # simplest starting value--but should in the end use another
  p.new <- pooledbinom.mle(x,m,n) #sum(x)/sum(m*n)
  #p.new <- pooledbinom.mle(x,m,n)
  done <- FALSE
  iter <- 0
  while(!done){
    iter <- iter+1
    p.old <- max(0,p.new)
    pip <- a - r * (1-p.old)^m

    vi <- n*m^2*(a-pip)^2 / ((1-p.old)^2 * pip * (1-pip))
    wi <- vi / sum(vi)
    vi.prime <- -(n*m^2*(a-pip)^2*(2*(m-1)*pip*(1-pip) + m*(a-pip)*(1-2*pip)))/((1-p.old)^3*pip^2*(1-pip)^2)
    wi.prime <- (vi.prime * sum(vi) - vi * sum(vi.prime)) / sum(vi)^2

    tSs <- sum(m*(a-pip)*(x-n*pip)/(pip*(1-pip)) )
    tSsp <- sum((m^2/(1-p.old)) *
                  ((x-n*(1-a))*pip^3 + (n*a*(1-a)-3*a*x)*pip^2 + a*(2*a+1)*x*pip - a^2*x) / (pip^2*(1-pip)^2) )
    p.new <- max(1e-8,p.old - tSs / tSsp) # 0 rather than 0.01 gave problems in some cases
    if(iter>max.iter){
      stop("Too many iterations")
    }
    if(abs(p.new-p.old)<tol)
      done <- 1
  }
  p.new
}


# Newton-Raphson iteration - so don't need root-finder
"ipIR.firth" <- "ipooledbinom.firth" <- function(x,m,n=rep(1,length(m)),
                                    sens=rep(1,length(m)),spec=rep(1,length(m)),
                                    tol=1e-12, max.iter=10000, p.start=NULL){
  # a, b match the manuscript; but sens, spec are more intuitive as function arguments
  a <- sens
  b <- spec
  # if really just the regular binomial, return the simple proportion
  if(all(m==1) & all(a == 1) & all(b == 1)) return(sum(x)/sum(n))
  if(sum(x)==0) return(0)

  r <- a + b - 1

  # Compute a (default) starting value:
  # When a starting value is not specified with p.start, this code calls this same function using
  # the Recall() functionality, but with forced a = 1 and b = 1.
  # This approach is to avoid having to reference a different function that gives the the perfect Firth
  # estimate, noting that this function itself gives the correct, perfect-test Firth estimate
  #  if a = 1 and b = 1.
  # So:  use p.start if it is explicitly specified, else compute a starting value
  if(!is.null(p.start)){
    p.new <- p.start
  } else {
    p.new <- NULL
    if(is.null(p.start) & is.null(p.new)){
      # if all pools are positive, the default starting value given below is 0, so make it something
      # positive
      if(sum(x) == sum(n)){
        p.new <- 1/sum(m*n)
      } else {
        N <- sum(n * m)
        mmw <- N / sum(n) # average pool size, N/sum(n)
        # uses proportion of positive pools, sum(x)/sum(n)
        p.new <- Recall(x,m,n,sens=1,spec=1,tol=tol,p.start=1 - (1 - sum(x)/sum(n))^(1/mmw)) #sum(x)/sum(m*n)
      }
    }
  }

  # now do the iterations
  done <- FALSE
  iter <- 0
  while(!done){
    iter <- iter+1
    #p.old <- max(0,p.new)
    p.old <- p.new
    pip <- a - r * (1-p.old)^m
    vi <- n*m^2*(a-pip)^2 / ((1-p.old)^2 * pip * (1-pip))
    wi <- vi / sum(vi)
    vi.prime <- -(n*m^2*(a-pip)^2*(2*(m-1)*pip*(1-pip) + m*(a-pip)*(1-2*pip)))/((1-p.old)^3*pip^2*(1-pip)^2)
    wi.prime <- (vi.prime * sum(vi) - vi * sum(vi.prime)) / sum(vi)^2
    tSs <- sum(m*(a-pip)*(x-n*pip)/(pip*(1-pip)) - 0.5*wi*(m-1))
    tSsp <- sum((m^2/(1-p.old)) *
                  ((x-n*(1-a))*pip^3 + (n*a*(1-a)-3*a*x)*pip^2 + a*(2*a+1)*x*pip - a^2*x) / (pip^2*(1-pip)^2) - 0.5*wi.prime*(m-1))
    #tSsp <- sum((m^2/(1)) *
    #              ((x-n*(a-pip))*pip^3 + (n*a*(1-a)-3*a*x)*pip^2 + a*(2*a+1)*x*pip - a^2*x) / (pip^2*(1-pip)^2) - 0.5*wi.prime*(m-1))/ (1-p.old)
    p.new <- max(0, p.old - tSs / tSsp) # 0 rather than 0.01 gave problems in some cases
    if(iter > max.iter){
      stop("Too many iterations")
    }
    if(abs(p.new-p.old)<tol)
      done <- 1
  }
  p.new
}


