
#################################################################################
# two-sample functions
#################################################################################
"pooledbinom.diff.loglike" <-
  function(d,s,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
    pooledbinom.loglike((s+d)/2,x1,m1,n1) + pooledbinom.loglike((s-d)/2,x2,m2,n2)


"score.p" <-
  function(p, x, m, n = rep(1,length(m)))
  {
    if(sum(x) == 0 & p > 0) return(-sum(m*n)/(1-p))
    if(p < 0 | p >= 1) return(0)
    sum(m*x/(1-(1-p)^m) - m*n)/(1-p)
  }	


"score2" <-
  function(s,d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
  {
    0.05*(score.p((s+d)/2,x1,m1,n1) + score.p((s-d)/2,x2,m2,n2))
  }

"score2.vec" <- 
  function(s,d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
  {
    ns <- length(s)
    s2 <- vector(length=ns)
    for(i in 1:ns)
      s2[i] <- score2(s[i],d,x1,m1,x2,m2,n1,n2)
    s2
  }

# VITAL here to have tol = something very small!!!!
"Shatd" <- 
  function(d, x1, m1, x2, m2, n1, n2)
  {
    #if(sum(x1)==0 & sum(x2)==0){
    #	N1 <- sum(m1*n1)
    #	N2 <- sum(m2*n2)
    #	return(2 - abs(N1-N2)*abs(d)/(N1+N2))
    #}
    if(sum(x1) == 0){
      if(d < 0){
        # g is value of dScore(s,d)/ds on the edge s = -d
        # must be vectorized for uniroot()
        g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
        {
          nd <- length(d)
          gv <- vector(length=nd)
          for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
          gv
        }
        brkt <- bracket.bounded.root(g,lower=-0.5,lbnd=-1,upper=-0.001,ubnd=0,
                                     x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
        dstar <- uniroot(g,lower = brkt[1],upper = brkt[2], tol = .Machine$double.eps^0.75,
                         x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        #			dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
        #							tol = .Machine$double.eps^0.75,
        #						x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        if(d <= dstar){ 
          return(-d)
        } else {
          
          brkt <- bracket.bounded.root(score2,lower= -0.5,lbnd=-d,upper=0,ubnd=2+d,
                                       d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
          ans <- uniroot(score2,lower = brkt[1],upper = brkt[2],
                         tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
          
          #				ans <- uniroot(score2,lower = -d + .Machine$double.eps^0.5,upper = 2 + d - .Machine$double.eps^0.5,
          #						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
          return(ans)
        }
      } 
      else {
        
        brkt <- bracket.bounded.root(score2,lower= d+0.001,lbnd=d,upper=2-d-0.001,ubnd=2-d,
                                     d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
        ans <- uniroot(score2,lower= brkt[1],upper = brkt[2],
                       tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        
        #			ans <- uniroot(score2,lower= d + .Machine$double.eps^0.5,upper = 2 - d -.Machine$double.eps^0.5,
        #						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        return(ans)
      }
    }
    if(sum(x2) == 0){
      if(d > 0){
        # g is value of dScore(s,d)/ds on the edge s = d
        # must be vectorized for uniroot()
        g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
        {
          nd <- length(d)
          gv <- vector(length=nd)
          for(i in 1:nd) gv[i] <- score2(d[i],d[i],x1,m1,x2,m2,n1,n2)
          gv
        }
        brkt <- bracket.bounded.root(g,lower=0.001,lbnd=0,upper=0.5,ubnd=1,
                                     x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
        dstar <- uniroot(g,lower = brkt[1],upper = brkt[2],
                         tol = .Machine$double.eps^0.75,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        #				dstar <- uniroot(g,lower = 0 + .Machine$double.eps^0.5,upper = 1 -.Machine$double.eps^0.5,
        #						tol = .Machine$double.eps^0.75,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        if(d >= dstar) 
          return(d)
        else {
          
          brkt <- bracket.bounded.root(score2,lower=d+0.001,lbnd=d,upper=2-d-0.001,ubnd=2-d,
                                       d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
          ans <- uniroot(score2,lower=brkt[1],upper = brkt[2],
                         tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
          
          #				ans <- uniroot(score2,lower= d + .Machine$double.eps,upper = 2 - d - .Machine$double.eps,
          #						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
          return(ans)
        }
      } 
      else {
        brkt <- bracket.bounded.root(score2,lower=-d+0.001,lbnd=-d,upper=2+d-0.001,ubnd=2+d,
                                     d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
        ans <- uniroot(score2,lower= brkt[1],upper=brkt[2],
                       tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        #			ans <- uniroot(score2,lower= -d+.Machine$double.eps^0.5,upper=2+d-.Machine$double.eps^0.5,
        #						tol = .Machine$double.eps^0.75,d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        return(ans)
      }
    }
    if(d < 0){
      brkt <- bracket.bounded.root(score2,lower=-d+0.001,lbnd=-d,upper=2+d-0.001,ubnd=2+d,
                                   d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
      ans <- uniroot(score2, lower =  brkt[1],
                     upper = brkt[2], 
                     tol = .Machine$double.eps^0.75, d = d, x1 = x1,
                     m1 = m1, x2 = x2, m2 = m2, n1 = n1,
                     n2 = n2)$root	
      
      #		ans <- uniroot(score2, lower =  - d + .Machine$double.eps^0.5,
      #					upper = 2 + d - .Machine$double.eps^0.5, 
      #					tol = .Machine$double.eps^0.75, d = d, x1 = x1,
      #					m1 = m1, x2 = x2, m2 = m2, n1 = n1,
      #					n2 = n2)$root
    }
    else {
      brkt <- bracket.bounded.root(score2,lower=d+0.025,lbnd=d,upper=2-d-0.05,ubnd=2-d,
                                   d=d,x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)
      ans <- uniroot(score2, lower = brkt[1], upper = brkt[2],
                     tol = .Machine$double.eps^0.75,d = d, x1 = x1,	m1 = m1, x2 = x2, m2 = m2, n1 = n1,	n2 = n2)$root
      #		ans <- uniroot(score2, lower = d + .Machine$double.eps^0.5, upper = 2 - d - .Machine$double.eps^0.5,
      #					tol = .Machine$double.eps^0.75,d = d, x1 = x1,	m1 = m1, x2 = x2, m2 = m2, n1 = n1,	n2 = n2)$root
    }
    ans
  }

"Shatd.vec" <- 
  function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
  {
    nd <- length(d)
    shat <- vector(length=nd)
    dhat <- pooledbinom.diff.mle(x1,m1,x2,m2,n1,n2)
    for(i in 1:nd){
      shat[i] <- Shatd(d[i],x1,m1,x2,m2,n1,n2)
    }
    shat
  }

# This is the profile likelihood ratio interval
"pooledbinom.diff.lrt.ci" <-
  function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha=0.05)
  {
    p1 <- pooledbinom.mle(x1, m1, n1)
    p2 <- pooledbinom.mle(x2, m2, n2)
    dhat <- p1 - p2
    shat <- p1 + p2
    
    if(sum(x1)==0 & sum(x2) == 0){
      # This is from Newcombe, p. 889
      # NOTE: exp(-3.84/2) = 0.1465
      N1 <- sum(m1*n1)
      N2 <- sum(m2*n2)
      g <- exp(-qchisq(1-alpha,1)/2) 
      lcl <- -1 + g^(1/N1)
      ucl <- 1 - g^(1/N2)
      return(c(p1=p1,p2=p2,diff=dhat,lcl=lcl,ucl=ucl,alpha=alpha))
    }
    
    loglike.p <- function(p,x,m,n){
      if(p > 0) sum(x*log((1-(1-p)^m))) + log(1-p)*sum(m*(n-x))
      else 0
    }
    
    loglike.ds <- function(d,s,x1,m1,x2,m2,n1,n2)
      loglike.p((s+d)/2,x1,m1,n1) + loglike.p((s-d)/2,x2,m2,n2)
    
    lrt.stat <- function(d0,dhat,shat,x1,m1,x2,m2,n1,n2)
    {
      shat0 <- Shatd(d0,x1,m1,x2,m2,n1,n2)
      2*(loglike.ds(dhat,shat,x1,m1,x2,m2,n1,n2) -
           loglike.ds(d0,shat0,x1,m1,x2,m2,n1,n2)) 
    }
    
    f <- function(d0,dhat,shat,x1,m1,x2,m2,n1,n2,alpha)
      lrt.stat(d0,dhat,shat,x1,m1,x2,m2,n1,n2) - qchisq(1-alpha,1)
    
    f.dhat <- f(dhat,dhat,shat,x1,m1,x2,m2,n1,n2,alpha)	 
    stepsize <- 10
    for(i in 1:stepsize){
      if(dhat - i/stepsize <= -1){
        lower.lim <- -1 + .Machine$double.eps
        break
      }
      if(f(dhat - i/stepsize,dhat,shat,x1,m1,x2,m2,n1,n2,alpha)*f.dhat < 0){
        lower.lim <- dhat - i/stepsize
        break
      }
    }
    # should have upper = dhat.mpl 
    
    lower <- uniroot(f, lower = lower.lim, upper = dhat , dhat=dhat,shat=shat,x1 = x1, m1 = 
                       m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,alpha=alpha)$root
    
    for(i in 1:stepsize){
      if(dhat + i/stepsize >= 1){
        upper.lim <- 1 - .Machine$double.eps
        break
      }
      if(f(dhat + i/stepsize,dhat,shat,x1,m1,x2,m2,n1,n2,alpha)*f.dhat < 0){
        upper.lim <- dhat + i/stepsize
        break
      }
    }
    # should have lower = dhat.mpl 
    upper <- uniroot(f, lower = dhat, upper = upper.lim, dhat=dhat,shat=shat,x1 = x1, m1 = 
                       m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,alpha=alpha)$root				
    
    c(p1=p1,p2=p2,diff=dhat,lower=lower,upper=upper,alpha=alpha)
  }

"pooledbinom.diff.lrtstat" <-
  function(d,x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha=0.05)
  {
    p1 <- pooledbinom.mle(x1, m1, n1)
    p2 <- pooledbinom.mle(x2, m2, n2)
    dhat <- p1 - p2
    shat <- p1 + p2
    
    loglike.p <- function(p,x,m,n){
      if(p > 0) sum(x*log((1-(1-p)^m))) + log(1-p)*sum(m*(n-x))
      else 0
    }
    
    loglike.ds <- function(d,s,x1,m1,x2,m2,n1,n2)
      loglike.p((s+d)/2,x1,m1,n1) + loglike.p((s-d)/2,x2,m2,n2)
    
    
    lrt.stat <- function(d0,dhat,shat,x1,m1,x2,m2,n1,n2)
    {
      shat0 <- Shatd(d0,x1,m1,x2,m2,n1,n2)
      2*(loglike.ds(dhat,shat,x1,m1,x2,m2,n1,n2) -
           loglike.ds(d0,shat0,x1,m1,x2,m2,n1,n2)) 
    }
    lrt.stat(d,dhat,shat,x1,m1,x2,m2,n1,n2)
  }

"pooledbinom.diff.lrtstat.vec" <-
  function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
  {
    n <- length(d)
    lmp <- vector(length=n)
    for(i in 1:n){
      lmp[i] <- pooledbinom.diff.lrtstat(d[i],x1,m1,x2,m2,n1,n2)
    }
    lmp
  }


# maximum likelihood estimate
"pooledbinom.diff.mle" <-
  function(x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))	
    pooledbinom.mle(x1,m1,n1) - pooledbinom.mle(x2,m2,n2)

"pooledbinom.diff.newcombe.ci" <-
  function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha = 0.05)
  {
    s1 <- as.vector(pooledbinom.score.ci(x1, m1, n1, alpha = alpha))
    # as.vector to drop names
    p1 <- s1[1]
    L1 <- s1[2]
    U1 <- s1[3]
    s2 <- as.vector(pooledbinom.score.ci(x2, m2, n2, alpha = alpha))
    p2 <- s2[1]
    L2 <- s2[2]
    U2 <- s2[3]
    thetahat <- p1 - p2
    delta <- sqrt((p1 - L1)^2 + (U2 - p2)^2)
    epsilon <- sqrt((U1 - p1)^2 + (p2 - L2)^2)
    c(p1 = p1, p2 = p2, diff = thetahat, lower = thetahat - delta, 
      upper = thetahat + epsilon, alpha = alpha)
  }

"pooledbinom.diff.newcombecscore.ci" <-
  function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha = 0.05)
  {
    s1 <- as.vector(pooledbinom.cscore.ci(x1, m1, n1, alpha = alpha))
    # as.vector to drop names
    p1 <- s1[1]
    L1 <- s1[2]
    U1 <- s1[3]
    s2 <- as.vector(pooledbinom.cscore.ci(x2, m2, n2, alpha = alpha))
    p2 <- s2[1]
    L2 <- s2[2]
    U2 <- s2[3]
    thetahat <- p1 - p2
    delta <- sqrt((p1 - L1)^2 + (U2 - p2)^2)
    epsilon <- sqrt((U1 - p1)^2 + (p2 - L2)^2)
    c(p1 = p1, p2 = p2, diff = thetahat, lower = thetahat - delta, 
      upper = thetahat + epsilon, alpha = alpha)
  }

"pooledbinom.diff.wald.ci" <-
  function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha = 0.05)
  {
    p1 <- pooledbinom.mle(x1, m1, n1)
    p1.var <- pooledbinom.mle.var(p1, m1, n1)
    p2 <- pooledbinom.mle(x2, m2, n2)
    p2.var <- pooledbinom.mle.var(p2, m2, n2)
    diff <- p1 - p2
    diff.var <- p1.var + p2.var
    z <- qnorm(1 - alpha/2)
    c(p1 = p1, p2 = p2, diff = diff, lower = diff - z * sqrt(diff.var),
      upper = diff + z * sqrt(diff.var), alpha = alpha)
  }


#######################################################################################################################
#######################################################################################################################
#######################################################
# Score methods
#######################################################
#######################################################################################################################
#######################################################################################################################
Z <- function(d, s, x1, m1, x2, m2, n1=rep(1,length(x1)), n2=rep(1,length(x2)))
{
  0.5*(score.p((d + s)/2, x1, m1, n1) - score.p((s -
                                                   d)/2, x2, m2, n2)) * sqrt((pooledbinom.mle.var(
                                                     (d + s)/2, m1, n1) + pooledbinom.mle.var((
                                                       s - d)/2, m2, n2)))
}

"score.p" <-
  function(p, x, m, n = rep(1,length(m)))
  {
    if(sum(x) == 0 & p >= 0) return(-sum(m*n)/(1-p))
    if(p <= 0 | p >= 1) return(0)
    sum(m*x/(1-(1-p)^m) - m*n)/(1-p)
  }

"score2" <-
  function(s,d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
  {
    0.5*(score.p((s+d)/2,x1,m1,n1) + score.p((s-d)/2,x2,m2,n2))
  }

"score2.vec" <-
  function(s,d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
  {
    ns <- length(s)
    s2 <- vector(length=ns)
    for(i in 1:ns) s2[i] <- score2(s[i],d,x1,m1,x2,m2,n1,n2)
    s2
  }

"pooledbinom.diff.scorestat" <-
  function(d, x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)))
  {
    dhat <- pooledbinom.diff.mle(x1,m1,x2,m2,n1,n2)
    
    Z <- function(d, shat, x1, m1, x2, m2, n1, n2)
    {
      if(sum(x1) == 0){
        g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
        {
          nd <- length(d)
          gv <- vector(length=nd)
          for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
          gv
        }
        dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
                         tol = .Machine$double.eps,
                         x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        N1 <- sum(m1*n1)
        
        if(d <= dstar){
          0.5*(-sum(m1*n1)/(1-(shat+d)/2) * as.numeric(shat + d > 0)  - score.p((shat - d)/2, x2, m2, n2)) * 
            sqrt(((1-(shat+d)/2)^2 / N1)  * as.numeric(shat+d > 0) + 
                   pooledbinom.mle.var((shat - d)/2, m2, n2))
          
        } else {
          0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
            sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2))
        }
      } else {
        #
        # This is the REGULAR (easy) case -- x <> 0
        #
        p1 <- (shat+d)/2
        p2 <- (shat-d)/2
        I1 <- pooledbinom.I(p1,m1,n1)
        I2 <- pooledbinom.I(p2,m2,n2)
        R <- I1 / (I1 + I2)
        V1 <- 1/I1
        V2 <- 1/I2
        S1 <- score.p(p1,x1,m1,n1)
        S2 <- score.p(p2,x2,m2,n2)
        ans <- (R*S1 - (1-R)*S2) * sqrt(V1 + V2)
        #cprint(ans)
        return(ans)
        #	0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
        #	sqrt((pooledbinom.mle.var((shat + d)/2, m1, n1) + pooledbinom.mle.var((
        #	shat - d)/2, m2, n2)))
      }
    }
    
    shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
    Z(d, shat, x1, m1, x2, m2, n1, n2)  # - qnorm(1-0.05/2)
  }


"pooledbinom.diff.scorestat.vec" <-
  function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
  {
    n <- length(d)
    sc <- vector(length=n)
    for(i in 1:n) {
      sc[i] <- pooledbinom.diff.scorestat(d[i],x1,m1,x2,m2,n1,n2)
    }
    sc
  }

"pooledbinom.diff.score.ci" <-
  function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha=0.05)
  {
    p1 <- pooledbinom.mle(x1, m1, n1)
    p2 <- pooledbinom.mle(x2, m2, n2)
    dhat <- p1 - p2
    shat <- p1 + p2
    
    #if(p1 == 0 & p2 == 0){
    if(sum(x1) == 0 & sum(x2)==0){
      # Not quite sure of this...MISSING 1/2!?!
      N1 <- sum(m1*n1)
      N2 <- sum(m2*n2)
      z2 <- qnorm(1-alpha/2)^2
      #lower <- -z2/(N2 + z2)
      lower <- -as.vector(pooledbinom.score.ci(x2,m2,n2,alpha=alpha)[3]) # -upper limit from pop'n 2
      #upper <- z2/(N1 + z2)
      upper <- as.vector(pooledbinom.score.ci(x1,m1,n1,alpha=alpha)[3]) # upper limit from pop'n 1
      return(c(p1=0,p2=0,diff=0,lower=lower,upper=upper,alpha=alpha))
    }
    if(sum(x1) == sum(n1) | sum(x2) == sum(n2))
      return(c(p1=pooledbinom.mle(x1,m1,n1),p2=pooledbinom.mle(x2,m2,n2),diff=0,lower=-1,upper=1))
    
    Z <- function(d, shat, x1, m1, x2, m2, n1, n2)
    {
      if(sum(x1) == 0){
        g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
        {
          nd <- length(d)
          gv <- vector(length=nd)
          for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
          gv
        }
        dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
                         tol = .Machine$double.eps,
                         x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        N1 <- sum(m1*n1)
        if(d <= dstar){
          return(0.5*(score.p((shat+d)/2, x1, m1, n1) * as.numeric(shat + d > 0)  - score.p((shat - d)/2, x2, m2, n2)) * 
                   sqrt(((1-(shat+d)/2)^2 / N1)  * as.numeric(shat+d > 0) + 
                          pooledbinom.mle.var((shat - d)/2, m2, n2)))
          
        } else {
          return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
                   sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
        }
      }
      if(sum(x2) == 0){
        g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
        {
          nd <- length(d)
          gv <- vector(length=nd)
          for(i in 1:nd) gv[i] <-  score2(d[i],d[i],x1,m1,x2,m2,n1,n2)
          gv
        }
        dstar <- uniroot(g,lower = 0 + .Machine$double.eps^0.5,upper = 1 -.Machine$double.eps^0.5,
                         tol = .Machine$double.eps,
                         x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        N2 <- sum(m2*n2)
        if(d >= dstar){
          return(0.5*(score.p((shat+d)/2, x1, m1, n1)  - 
                        score.p((shat - d)/2, x2, m2, n2)*as.numeric(shat - d > 0)) * 
                   sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + (1-(shat-d)/2)^2/N2 * as.numeric(shat-d>0)))
          
        } else {
          return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
                   sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
        }
      }
      #
      # This is the REGULAR (easy) case -- x <> 0
      #
      0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
        sqrt((pooledbinom.mle.var((shat + d)/2, m1, n1) + pooledbinom.mle.var((
          shat - d)/2, m2, n2)))
    }
    
    f.lower <- function(d, x1, m1, x2, m2, n1, n2,dhat)
    {
      shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
      z <- Z(d,shat,x1,m1,x2,m2,n1,n2) - qnorm(1-alpha/2)
      z
    }
    f.upper <- function(d, x1, m1, x2, m2, n1, n2,dhat)
    {
      shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
      z <- Z(d,shat,x1,m1,x2,m2,n1,n2) + qnorm(1-alpha/2)
      z
    }
    
    # Bound the root, in a rather slow but safe way
    f.dhat <- f.lower(dhat,x1,m1,x2,m2,n1,n2,dhat)	 
    
    stepsize <- 10
    for(i in 1:stepsize){
      if(dhat - i/stepsize <= -1){
        lower.lim <- -1 + .Machine$double.eps
        break
      }
      if(f.lower(dhat - i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
        lower.lim <- dhat - i/stepsize
        break
      }
    }
    lower <- uniroot(f.lower, lower = lower.lim, upper = dhat , x1 = x1, m1 = 
                       m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root
    f.dhat <- f.upper(dhat,x1,m1,x2,m2,n1,n2,dhat)	 
    for(i in 1:stepsize){
      if(dhat + i/stepsize >= 1){
        upper.lim <- 1 - .Machine$double.eps
        break
      }
      if(f.upper(dhat + i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
        upper.lim <- dhat + i/stepsize
        break
      }
    }
    #if(!exists("upper.lim")) upper.lim <- 1-.Machine$double.eps^0.5
    
    upper <- uniroot(f.upper, lower = dhat, upper = upper.lim,x1 = x1, m1 = 
                       m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root
    c(p1 = p1, p2 = p2, diff = dhat, lower = lower, upper = upper, alpha=alpha)
  }


"pooledbinom.diff.cscore.ci" <-
  function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha=0.05)
  {
    p1 <- pooledbinom.mle(x1, m1, n1)
    p2 <- pooledbinom.mle(x2, m2, n2)
    dhat <- p1 - p2
    shat <- p1 + p2
    chi2 <- qnorm(1 - alpha/2)^2
    if(sum(x1) == 0 & sum(x2)==0){
      # Not quite sure of this...MISSING 1/2!?!
      N1 <- sum(m1*n1)
      N2 <- sum(m2*n2)
      z2 <- qnorm(1-alpha/2)^2
      lower <- -z2/(N2 + z2)
      upper <- z2/(N1 + z2)
      return(c(p1=0,p2=0,diff=0,lower=lower,upper=upper,alpha=alpha))
    }
    if(sum(x1) == sum(n1) | sum(x2) == sum(n2))
      return(c(p1=pooledbinom.mle(x1,m1,n1),p2=pooledbinom.mle(x2,m2,n2),diff=0,lower=-1,upper=1))
    
    Z <- function(d, shat, x1, m1, x2, m2, n1, n2)
    {
      if(sum(x1) == 0){
        g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
        {
          nd <- length(d)
          gv <- vector(length=nd)
          for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
          gv
        }
        dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
                         tol = .Machine$double.eps,
                         x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        N1 <- sum(m1*n1)
        if(d <= dstar){
          return(0.5*(score.p((shat+d)/2, x1, m1, n1) * as.numeric(shat + d > 0)  - score.p((shat - d)/2, x2, m2, n2)) * 
                   sqrt(((1-(shat+d)/2)^2 / N1)  * as.numeric(shat+d > 0) + 
                          pooledbinom.mle.var((shat - d)/2, m2, n2)))
          
        } else {
          return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
                   sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
        }
      }
      if(sum(x2) == 0){
        g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
        {
          nd <- length(d)
          gv <- vector(length=nd)
          for(i in 1:nd) gv[i] <-  score2(d[i],d[i],x1,m1,x2,m2,n1,n2)
          gv
        }
        dstar <- uniroot(g,lower = 0 + .Machine$double.eps^0.5,upper = 1 -.Machine$double.eps^0.5,
                         tol = .Machine$double.eps,
                         x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        N2 <- sum(m2*n2)
        if(d >= dstar){
          return(0.5*(score.p((shat+d)/2, x1, m1, n1)  - 
                        score.p((shat - d)/2, x2, m2, n2)*as.numeric(shat - d > 0)) * 
                   sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + (1-(shat-d)/2)^2/N2 * as.numeric(shat-d>0)))
          
        } else {
          return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
                   sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
        }
      }
      #
      # This is the REGULAR (easy) case -- x <> 0
      #
      0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
        sqrt((pooledbinom.mle.var((shat + d)/2, m1, n1) + pooledbinom.mle.var((
          shat - d)/2, m2, n2)))
    }
    
    mu3 <- function(d,shat,m1,m2,n1,n2)
    {
      I1 <- pooledbinom.I((shat+d)/2,m1,n1)
      I2 <- pooledbinom.I((shat-d)/2,m2,n2)
      R <- I1/(I1+I2)
      (1-R)^3*pooledbinom.mu3((shat+d)/2,m1,n1) - R^3*pooledbinom.mu3((shat-d)/2,m2,n2)
    }
    gamma1 <- function(d,shat,m1,m2,n1,n2)
    {
      mu3(d,shat,m1,m2,n1,n2)*
        (pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat-d)/2,m2,n2))^(3/2)
    }
    
    f.lower <- function(d, x1, m1, x2, m2, n1, n2,dhat)
    {
      shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
      z <- Z(d,shat,x1,m1,x2,m2,n1,n2) - gamma1(d,shat,m1,m2,n1,n2)*(chi2 - 1)/6 - qnorm(1-alpha/2)
      z
    }
    f.upper <- function(d, x1, m1, x2, m2, n1, n2,dhat)
    {
      shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
      z <- Z(d,shat,x1,m1,x2,m2,n1,n2) - gamma1(d,shat,m1,m2,n1,n2)*(chi2 - 1)/6 + qnorm(1-alpha/2)
      z
    }
    
    # Bound the root, in a rather slow but safe way
    f.dhat <- f.lower(dhat,x1,m1,x2,m2,n1,n2,dhat)	 
    
    stepsize <- 10
    for(i in 1:stepsize){
      if(dhat - i/stepsize <= -1){
        lower.lim <- -1 + .Machine$double.eps
        break
      }
      if(f.lower(dhat - i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
        lower.lim <- dhat - i/stepsize
        break
      }
    }
    lower <- uniroot(f.lower, lower = lower.lim, upper = dhat , x1 = x1, m1 = 
                       m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root
    
    f.dhat <- f.upper(dhat,x1,m1,x2,m2,n1,n2,dhat)	 
    for(i in 1:stepsize){
      if(dhat + i/stepsize >= 1){
        upper.lim <- 1 - .Machine$double.eps
        break
      }
      if(f.upper(dhat + i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
        upper.lim <- dhat + i/stepsize
        break
      }
    }
    
    upper <- uniroot(f.upper, lower = dhat, upper = upper.lim,x1 = x1, m1 = 
                       m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root  		 
    
    
    c(p1=p1,p2=p2,diff = dhat, lower = lower, upper = upper,alpha=alpha)
  }


"pooledbinom.diff.bcscore.ci" <-
  function(x1, m1, x2, m2, n1 = rep(1, length(x1)), n2 = rep(1, length(x2)),alpha=0.05)
  {
    p1 <- pooledbinom.mle(x1, m1, n1)
    p2 <- pooledbinom.mle(x2, m2, n2)
    dhat <- p1 - p2
    shat <- p1 + p2
    chi2 <- qnorm(1 - alpha/2)^2
    if(sum(x1) == 0 & sum(x2)==0){
      # Not quite sure of this...MISSING 1/2!?!
      N1 <- sum(m1*n1)
      N2 <- sum(m2*n2)
      z2 <- qnorm(1-alpha/2)^2
      lower <- -z2/(N2 + z2)
      upper <- z2/(N1 + z2)
      return(c(p1=0,p2=0,diff=0,lower=lower,upper=upper,alpha=alpha))
    }
    if(sum(x1) == sum(n1) | sum(x2) == sum(n2))
      return(c(p1=pooledbinom.mle(x1,m1,n1),p2=pooledbinom.mle(x2,m2,n2),diff=0,lower=-1,upper=1))
    
    Z <- function(d, shat, x1, m1, x2, m2, n1, n2)
    {
      if(sum(x1) == 0){
        g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
        {
          nd <- length(d)
          gv <- vector(length=nd)
          for(i in 1:nd) gv[i] <-  score2(-d[i],d[i],x1,m1,x2,m2,n1,n2)
          gv
        }
        dstar <- uniroot(g,lower = -1 + .Machine$double.eps^0.5,upper = 0 -.Machine$double.eps^0.5,
                         tol = .Machine$double.eps,
                         x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        N1 <- sum(m1*n1)
        if(d <= dstar){
          return(0.5*(score.p((shat+d)/2, x1, m1, n1) * as.numeric(shat + d > 0)  - score.p((shat - d)/2, x2, m2, n2)) * 
                   sqrt(((1-(shat+d)/2)^2 / N1)  * as.numeric(shat+d > 0) + 
                          pooledbinom.mle.var((shat - d)/2, m2, n2)))
          
        } else {
          return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
                   sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
        }
      }
      if(sum(x2) == 0){
        g <- function(d,x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)))
        {
          nd <- length(d)
          gv <- vector(length=nd)
          for(i in 1:nd) gv[i] <-  score2(d[i],d[i],x1,m1,x2,m2,n1,n2)
          gv
        }
        dstar <- uniroot(g,lower = 0 + .Machine$double.eps^0.5,upper = 1 -.Machine$double.eps^0.5,
                         tol = .Machine$double.eps,
                         x1=x1,m1=m1,x2=x2,m2=m2,n1=n1,n2=n2)$root
        N2 <- sum(m2*n2)
        if(d >= dstar){
          return(0.5*(score.p((shat+d)/2, x1, m1, n1)  - 
                        score.p((shat - d)/2, x2, m2, n2)*as.numeric(shat - d > 0)) * 
                   sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + (1-(shat-d)/2)^2/N2 * as.numeric(shat-d>0)))
          
        } else {
          return(0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
                   sqrt(pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat - d)/2, m2, n2)))
        }
      }
      #
      # This is the REGULAR (easy) case -- x <> 0
      #
      0.5*(score.p((shat+d)/2, x1, m1, n1) - score.p((shat - d)/2, x2, m2, n2)) * 
        sqrt((pooledbinom.mle.var((shat + d)/2, m1, n1) + pooledbinom.mle.var((
          shat - d)/2, m2, n2)))
    }
    mu3 <- function(d,shat,m1,m2,n1,n2)
    {
      I1 <- pooledbinom.I((shat+d)/2,m1,n1)
      I2 <- pooledbinom.I((shat-d)/2,m2,n2)
      R <- I1/(I1+I2)
      (1-R)^3*pooledbinom.mu3((shat+d)/2,m1,n1) - R^3 * pooledbinom.mu3((shat-d)/2,m2,n2)
    }
    gamma1 <- function(d,shat,m1,m2,n1,n2)
    {
      mu3(d,shat,m1,m2,n1,n2)*
        (pooledbinom.mle.var((shat+d)/2,m1,n1) + pooledbinom.mle.var((shat-d)/2,m2,n2))^(3/2)
    }
    bias <- function(d,shat,m1,m2,n1,n2)
    {
      I1 <- pooledbinom.I((shat+d)/2,m1,n1)
      I2 <- pooledbinom.I((shat-d)/2,m2,n2)
      R <- I1/(I1+I2)
      R^(3/2) * sqrt(I2) * pooledbinom.bias((shat+d)/2,m1,n1) -
        (1-R)^(3/2) * sqrt(I1) * pooledbinom.bias((shat-d)/2,m2,n2)
    }
    
    f.lower <- function(d, x1, m1, x2, m2, n1, n2,dhat)
    {
      shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
      z <- Z(d,shat,x1,m1,x2,m2,n1,n2) - bias(d,shat,m1,m2,n1,n2) - gamma1(d,shat,m1,m2,n1,n2)*(chi2 - 1)/6 - qnorm(1-alpha/2)
      z
    }
    f.upper <- function(d, x1, m1, x2, m2, n1, n2,dhat)
    {
      shat <- Shatd(d, x1, m1, x2, m2, n1, n2)
      z <- Z(d,shat,x1,m1,x2,m2,n1,n2) - bias(d,shat,m1,m2,n1,n2) - gamma1(d,shat,m1,m2,n1,n2)*(chi2 - 1)/6 + qnorm(1-alpha/2)
      z
    }
    
    # Bound the root, in a rather slow but safe way
    f.dhat <- f.lower(dhat,x1,m1,x2,m2,n1,n2,dhat)	 
    
    stepsize <- 10
    for(i in 1:stepsize){
      if(dhat - i/stepsize <= -1){
        lower.lim <- -1 + .Machine$double.eps
        break
      }
      if(f.lower(dhat - i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
        lower.lim <- dhat - i/stepsize
        break
      }
    }
    lower <- uniroot(f.lower, lower = lower.lim, upper = dhat , x1 = x1, m1 = 
                       m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root
    
    f.dhat <- f.upper(dhat,x1,m1,x2,m2,n1,n2,dhat)	 
    for(i in 1:stepsize){
      if(dhat + i/stepsize >= 1){
        upper.lim <- 1 - .Machine$double.eps
        break
      }
      if(f.upper(dhat + i/stepsize,x1,m1,x2,m2,n1,n2,dhat)*f.dhat < 0){
        upper.lim <- dhat + i/stepsize
        break
      }
    }
    
    upper <- uniroot(f.upper, lower = dhat, upper = upper.lim,x1 = x1, m1 = 
                       m1, x2 = x2, m2 = m2, n1 = n1, n2 = n2,dhat=dhat)$root  		 
    
    c(p1=p1,p2=p2,diff = dhat, lower = lower, upper = upper,alpha=alpha)
  }


########################
# MIR
########################

"pooledbinom.diff.mir.ci" <-
  function(x1,m1,x2,m2,n1=rep(1,length(x1)),n2=rep(1,length(x2)),alpha=0.05)
  {
    N1 <- sum(m1*n1)
    N2 <- sum(m2*n2)
    mir1 <- sum(x1)/N1
    mir2 <- sum(x2)/N2
    mir1.var <- mir1*(1-mir1)/N1
    mir2.var <- mir2*(1-mir2)/N2
    mir.diff <- mir1 - mir2
    mir.diff.stderr <- sqrt(mir1.var + mir2.var)
    z <- qnorm(1-alpha/2)
    c(p1=mir1,p2=mir2,diff=mir.diff,lower=mir.diff - z * mir.diff.stderr,upper=mir.diff + z*mir.diff.stderr,alpha=alpha)
  }

