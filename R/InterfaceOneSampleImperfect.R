
"ipooledBin" <- function(x,...){
  UseMethod("ipooledBin")
}

"ipooledBin.default" <- function(x,m,n=rep(1,length(x)), group,
           pt.method = c("firth","mle"),
           sens=rep(1,length(x)),spec=rep(1,length(x)),  scale=1,tol=.Machine$double.eps^0.5,
           max.iter=10000,p.start=NULL,...) {
    call <- match.call()
    call[[1]] <- as.name("ipooledBin")
    pt.method <- match.arg(pt.method)
    if(!missing(group)){
      groups <- unique(group)
      nGroups <- length(groups)
      ans <- vector(mode="list",length=nGroups)

      for(i in 1:nGroups){
        switch(pt.method,
               firth = {
                ans[[i]] <- ipooledbinom.firth(x[group==groups[i]],
                                             m[group==groups[i]],
                                             n[group==groups[i]],
                                             sens=sens,spec=spec,tol=tol,
                                             max.iter=max.iter,p.start=p.start)
                },
                mle = {
                  ans[[i]] <- ipooledbinom.firth(x[group==groups[i]],
                                               m[group==groups[i]],
                                               n[group==groups[i]],
                                               sens=sens,spec=spec,tol=tol,
                                               max.iter=max.iter)
                  }
                  )

        ans[[i]] <- structure(ans[[i]],x=x[group==groups[i]],
                              m=m[group==groups[i]],
                              n=n[group==groups[i]],class="ipooledBIn",call=call)
        #class(ans[[i]]) <- "ipooledBin"
        #attr(ans[[i]],"call") <- call
      }
      #names(ans) <- groups
      #class(ans) <- "pooledBinList" # added this so it will be able to use the print.pooledBinList function""
      #attributes(ans) <- vars
      #class(ans) <- "pooledBinList"
      #attributes(ans) <- list(class = "pooledBinList", names = groups,
      #                        x.var = vars$x, m.var = vars$m, n.var = vars$n, group.var = vars$group)

      #names(ans) <- groups
      if(class(substitute(group)) == "name") group.var <- deparse(substitute(group))
      else group.var <- "Group"

      #attributes(ans) <- list(class = "pooledBinList", group.names = groups, group.var = group.var,call=call)

      ans <- structure(ans, class = "ipooledBinList", group.names = groups, group.var = group.var,
                       x=x,m=m,n=n,
                       sens=sens,spec=spec,scale = scale, pt.method=pt.method,call=call)

    } else {
      switch(pt.method,
             firth = {
      ans <- ipooledbinom.firth(x, m, n,
                              sens=sens,spec=spec,tol=tol,
                              max.iter=max.iter,p.start=p.start)
             },
      mle = {
      ans <- ipooledbinom.mle(x, m, n,
                              sens=sens,spec=spec,tol=tol,
                              max.iter=max.iter)

      }
      )
      #attributes(ans) <- list(attributes(ans), call=call)
      #names(ans) <- "1"
      #ans$call <- call
      #class(ans) <- "pooledBin"
      #structure(ans,class="pooledBin",call=call)
      #if(class(substitute(group)) == "name") group.var <- deparse(substitute(group))
      #else group.var <- "Group"
      ans <- structure(ans,class="ipooledBin",
                       x=x,m=m,n=n,
                       sens=sens,spec=spec,scale=scale, pt.method=pt.method,call=call,group.names = "", group.var = "")
    }
    ans
  }


"ipooledBin.formula" <- function(x, data,
           pt.method = c("firth","mle"),
           sens = 1, spec = 1,
           scale = 1, tol=.Machine$double.eps^0.5,
           max.iter=10000,p.start=NULL,...){
    call <- match.call()
    call[[1]] <- as.name("ipooledBin")

    pt.method <- match.arg(pt.method)

    if(missing(data))
      data <- environment(x)

    # omit records with missing data -- note, if data contains records missing
    # anywhere (even not in X, M, N, Group, they are omitted), so care should be
    # used in subsetting before the call to be asured of the desired analysis
    #data <- na.omit(data)

    vars <- pooledBinParseFormula(x, data)
    if(any(sapply(vars,length)>1)) stop("only variable names permitted in formula; perhaps use the default call")

    # retrieve values from data using the character name
    # use the eval(parse(text=XX), data) construct in case data is left unspecified,
    # and because we have the names of the variables available
    # from pooledBin.deparseFormula
    x <- eval(parse(text=vars$x), data)
    m <- eval(parse(text=vars$m), data)
    if(!is.null(vars$n))
      n <- eval(parse(text=vars$n), data)
    else n <- rep(1,length(x)) # default n

    if(!is.null(vars$group)){
      group <- eval(parse(text=vars$group), data)
      groups <- unique(group)
      nGroups <- length(groups)
    }  else{
      group <- rep("SINGLE",length(x))
      groups <- unique(group)
      nGroups <- 1
    }
    if(nGroups > 1){
      ans <- vector(mode="list",length=nGroups)
      for(i in 1:nGroups){
        switch(pt.method,
               firth = {
                 ans[[i]] <- ipooledbinom.firth(x[group==groups[i]],
                                              m[group==groups[i]],
                                              n[group==groups[i]],
                                              sens=sens,spec=spec,tol=tol,
                                              max.iter=max.iter,p.start=p.start)
               },
               mle = {
                 ans[[i]] <- ipooledbinom.mle(x[group==groups[i]],
                                              m[group==groups[i]],
                                              n[group==groups[i]],
                                              sens=sens,spec=spec,tol=tol,
                                              max.iter=max.iter)
               }
        )
        ans[[i]] <- structure(ans[[i]],x=x[group==groups[i]],
                              m=m[group==groups[i]],
                              n=n[group==groups[i]],class="ipooledBIn",call=call)
        #class(ans[[i]]) <- "ipooledBin"
        #attributes(ans[[i]]) <- list(class="pooledBin",call=call)
        #attr(ans[[i]],"call") <- call

      }
      names(ans) <- groups
      ans <- structure(ans, class = "ipooledBinList", group.names = groups,
                       x.var = vars$x, m.var = vars$m, n.var = vars$n, group.var = vars$group,
                       x=x,m=m,n=n,
                       scale = scale, sens=sens,spec=spec,pt.method = pt.method,call = call)
    } else {
      switch(pt.method,
             firth = {
               ans <- ipooledbinom.firth(x, m, n,
                                       sens=sens,spec=spec,tol=tol,
                                       max.iter=max.iter,p.start=p.start)
             },
             mle = {
               ans <- ipooledbinom.mle(x, m, n,
                                     sens=sens,spec=spec,tol=tol,
                                     max.iter=max.iter)

             }
      )
      #names(ans) <- "1"
      ans <- structure(ans, class = "ipooledBin", group.names = names(ans),
                       x.var = vars$x, m.var = vars$m, n.var = vars$n, group.var = vars$group,
                       x=x,m=m,n=n,
                       sens=sens,spec=spec,scale=scale,pt.method=pt.method,call = call)

    }
    #attributes(ans,"x.var") <- vars$x
    #attributes(ans,"m.var") <- vars$m
    #attributes(ans,"n.var") <- vars$n
    #if(!is.null(vars$group)) attributes(ans,"group.var") <- vars$group

    ans

  }



"print.ipooledBin" <- function(x, ..., scale=attr(x,"scale")){
    args <- list(...)
    if(is.null(scale)) scale <- 1
    if(is.null(args$digits)) digits <- 4
    else digits <- args$digits
    p <- round(scale*x,digits)
    mat <- matrix(c(p,scale),nrow=1) # really to match Hmisc's binconf()
    dimnames(mat) <- list(c(""),c("P","Scale"))
    if(scale == 1) mat <- mat[,-2,drop=FALSE]
    print(mat,...)
    invisible(x)
  }

"print.ipooledBinList" <- function(x, ...){
    n <- length(x)
    out <- data.frame(Group = attr(x,"group.names"),
                      PointEst = rep(0,n),
                      Scale = rep(1,n))
    if(!is.null(attr(x,"group.var"))) names(out)[1] <- attr(x,"group.var")
    for(i in 1:n){
      out[i,2:3] <-  attr(x,"scale") * c(x[[i]], 1)
    }
    if(all(out$Scale == 1)) out$Scale <- NULL # if scale = 1, don't bother printing (do print in summary, though)
    print(out,...)
    invisible(x)
  }


"as.data.frame.ipooledBin" <- function(x, row.names = NULL, optional = FALSE, ...){
  args <- list(...)
  if(is.null(scale)) scale <- 1
  if(is.null(args$digits)) digits <- 4
  else digits <- args$digits
  p <- round(scale*x,digits)
  mat <- matrix(c(p,scale),nrow=1) # really to match Hmisc's binconf()
  dimnames(mat) <- list(c(""),c("P","Scale"))
  if(scale == 1) mat <- mat[,-2,drop=FALSE]
  as.data.frame(mat)
}

"as.data.frame.ipooledBinList" <- function(x, row.names = NULL, optional = FALSE, ...){
  n <- length(x)
  out <- data.frame(Group = attr(x,"group.names"),
                    PointEst = rep(0,n),
                    Scale = rep(1,n))
  if(!is.null(attr(x,"group.var"))) names(out)[1] <- attr(x,"group.var")
  for(i in 1:n){
    out[i,2:3] <-  attr(x,"scale") * c(x[[i]], 1)
  }
  if(all(out$Scale == 1)) out$Scale <- NULL # if scale = 1, don't bother printing (do print in summary, though)
  out
}


"summary.ipooledBin" <-
  function(object, ...){
    x <- object
    args <- list(...)
    scale <- attr(x,"scale")
    structure(x,class="summary.ipooledBin",scale=scale,
                       x=attr(x,"x"),m=attr(x,"m"),n=attr(x,"n"),
              sens=attr(x,"sens"),spec=attr(x,"spec"),pt.method=attr(x,"pt.method"),call=attr(object,"call"))
  }

"print.summary.ipooledBin" <-
  function(x, ...){
    args <- list(...)
    scale <- attr(x,"scale")
    if(is.null(args$digits)) digits <- 4
    else digits <- args$digits
    cat("Estimation of Binomial Proportion for Pooled Data using an Imperfect Test\n\n")
    cat(paste0("Call: ", deparse(attr(x,"call"),width.cutoff = 120),"\n\n"))


    if(all(attr(x,"sens")==1)) sens <- 1
    else sens <- attr(x,"sens")
    if(all(attr(x,"spec")==1)) spec <- 1
    else spec <- attr(x,"spec")
    cat(paste0("Call: ", deparse(attr(x,"call"),width.cutoff = 120),"\n\n"))
    cat(paste0("Sensitivity : ",paste0(" ",sens,collapse=","),"\n"))
    cat(paste0("Specificity : ",paste0(" ",spec,collapse=","),"\n"))
    cat("\n")
    cat(paste("Point estimator:",attr(x,"pt.method"),"\n"))
    cat("\n")

    cat(paste("Number of individuals:",sum(attr(x,"n") * attr(x,"m")),"\n"))
    cat(paste("Number of pools:",sum(attr(x,"n")),"\n"))
    cat(paste("Number of positive pools:",sum(attr(x,"x")),"\n"))
    cat("\n")
    print.ipooledBin(x, ..., scale=scale)
    invisible(x)
  }


"summary.ipooledBinList" <- function(object, ...){
  grp.names <- as.character(attr(object,"group.names"))
  "sumf" <- function(x) c(P=x,N = sum(attr(x,"n") * attr(x,"m")), NumPools = sum(attr(x,"n")), NumPosPools = sum(attr(x,"x")))
  out <- lapply(object, sumf)
  #attributes(out) <- list(class = "summary.pooledBinList", names = attr(object,"names"),group.var = attr(object,"group.var"))
  #attributes(out) <- list(class = "summary.pooledBinList",
  #                        names = grp.names,group.var = attr(object,"group.var"))
  structure(out, class = "summary.ipooledBinList",
            names = grp.names,group.var = attr(object,"group.var"),
            scale = attr(object,"scale"),
                       x=attr(object,"x"),m=attr(object,"m"),n=attr(object,"n"),
            sens = attr(object,"sens"),
            spec = attr(object,"spec"),
            pt.method = attr(object,"pt.method"),
            call = attr(object,"call"))
}

"print.summary.ipooledBinList" <- function(x, ...){
  n <- length(x)
  out <- data.frame(Group = names(x),
                    PointEst = rep(0,n),
                    N = rep(0,n),
                    NumPools = rep(0,n),
                    NumPosPools = rep(0,n),
                    Scale = rep(1,n))
  if(!is.null(attr(x,"group.var"))) names(out)[1] <- attr(x,"group.var")
  for(i in 1:n)
    out[i,2:6] <- c(attr(x,"scale") * x[[i]][["P"]], x[[i]][["N"]], x[[i]][["NumPools"]], x[[i]][["NumPosPools"]], attr(x,"scale"))

  #if(is.null(digits)) digits <- 4
  #else digits <- args$digits
  cat("\n")
  cat("Estimation of Binomial Proportion for Pooled Data using an Imperfect Test\n")
  cat(paste0("Call: ", deparse(attr(x,"call"),width.cutoff = 120),"\n\n"))
  cat(paste0("\nCall: ", deparse(attr(x,"call"),width.cutoff=120),"\n\n"))
  if(all(attr(x,"sens")==1)) sens <- 1
  else sens <- attr(x,"sens")
  if(all(attr(x,"spec")==1)) spec <- 1
  else spec <- attr(x,"spec")
  cat(paste0("Sensitivity : ",paste0(" ",sens,collapse=","),"\n"))
  cat(paste0("Specificity : ",paste0(" ",spec,collapse=","),"\n"))
  cat(paste("Point estimator :",attr(x,"pt.method"),"\n"))

  print(out,...)
  invisible(x)
}


