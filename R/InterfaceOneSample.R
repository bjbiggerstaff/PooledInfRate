
"pooledBin" <- function(x,...){
  UseMethod("pooledBin")
}

"pooledBin.default" <- function(x,m,n=rep(1,length(x)), group,
           pt.method = c("firth","gart","bc-mle","mle","mir"),
           ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
           scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
    call <- match.call()
    call[[1]] <- as.name("pooledBin")

    xmn.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n))
    x <- x[xmn.nomiss]
    m <- m[xmn.nomiss]
    n <- n[xmn.nomiss]

    mn.nozero <- (m>0 & n>0)
    x <- x[mn.nozero]
    m <- m[mn.nozero]
    n <- n[mn.nozero]

    if(!missing(group)){
    xmng.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n) & !is.na(group))
    group <- group[xmng.nomiss]
      groups <- unique(group)
      nGroups <- length(groups)
      ans <- vector(mode="list",length=nGroups)

      for(i in 1:nGroups){
        ans[[i]] <- pooledBin.fit(x[group==groups[i]],
                                                    m[group==groups[i]],
                                                    n[group==groups[i]],
                                                    pt.method=pt.method,
                                                    ci.method=ci.method,
                                                    scale=scale,alpha=alpha,tol=tol)
        class(ans[[i]]) <- "pooledBin"
        attr(ans[[i]],"call") <- call
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

      #
      ans.lst <- ans
      ans <- data.frame(Group = attr(x, "group.names"),
                        PointEst = rep(0,  nGroups),
                        Lower = rep(0, nGroups),
                        Upper = rep(0, nGroups),
                        Scale = rep(1,  nGroups))
      if (!is.null(attr(ans.lst, "group.var")))
        names(ans)[1] <- attr(ans.lst, "group.var")
      for (i in 1:n) out[i, 2:5] <- ans.lst[[i]]$scale * c(ans.lst[[i]]$p,  ans.lst[[i]]$lcl, ans.lst[[i]]$ucl, 1)
      if (all(ans$Scale == 1)) ans$Scale <- NULL
print(ans)
      ans <- structure(ans, class = "pooledBinList", fullList = ans.lst, group.names = groups, group.var = group.var,scale=scale,call=call)

    } else {
      ans <- pooledBin.fit(x, m, n,
                           pt.method=pt.method,
                           ci.method=ci.method,
                           scale=scale,alpha=alpha,tol=tol)
      #attributes(ans) <- list(attributes(ans), call=call)
      #names(ans) <- "1"
      #ans$call <- call
      #class(ans) <- "pooledBin"
      #structure(ans,class="pooledBin",call=call)
      #if(class(substitute(group)) == "name") group.var <- deparse(substitute(group))
      #else group.var <- "Group"
      ans.lst <- ans
      ans <- data.frame(P=scale*ans.lst$p,
                        Lower = scale*ans.lst$lcl,
                        Upper = scale*ans.lst$ucl,
                        Scale = scale) # really to match Hmisc's binconf()
      if(scale == 1) ans$Scale <- NULL
      ans <- structure(ans,class="pooledBin",fullList = ans.lst,call=call,group.names = "", group.var = "", scale=scale)
    }
    ans
  }


"pooledBin.formula" <- function(x, data,
           pt.method = c("firth","gart","bc-mle","mle","mir"),
           ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
           scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...){
    call <- match.call()
    call[[1]] <- as.name("pooledBin")

    if(missing(data))
      data <- environment(x)

    vars <- pooledBinParseFormula(x, data)
    if(any(sapply(vars,length)>1)) stop("only variable names permitted in formula; perhaps use the default call")



    # omit records with missing data -- note, if data contains records missing
    # anywhere (even not in X, M, N, Group, they are omitted), so care should be
    # used in subsetting before the call to be assured of the desired analysis
    # restrict to only those variables needed before subsetting to avoid that issue.
    if(!missing(data)){
      # restrict the data to the variables needed before removing records with any NA
      # --this doesn't help when no data are specified, so that's dealt with below
      data <- data[,unlist(vars)]
      data <- na.omit(data)
    }


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

      # restrict to data with no missing x, m, n, group
    if(missing(data)){
      mn.nozero <- (m>0 & n>0)
      x <- x[mn.nozero]
      m <- m[mn.nozero]
      n <- n[mn.nozero]

      xmng.nomiss <- (!is.na(x) & !is.na(m) & !is.na(n) & !is.na(group))
      x <- x[xmng.nomiss]
      m <- m[xmng.nomiss]
      n <- n[xmng.nomiss]


      group <- group[xmng.nomiss]
      groups <- unique(group)
      nGroups <- length(groups)
    }

    mn.nozero <- (m>0 & n>0)
    x <- x[mn.nozero]
    m <- m[mn.nozero]
    n <- n[mn.nozero]
    group <- group[mn.nozero]
    groups <- unique(group)
    nGroups <- length(groups)

    if(nGroups > 1){
      ans <- vector(mode="list",length=nGroups)
      for(i in 1:nGroups){
        ans[[i]] <- pooledBin.fit(x[group==groups[i]],
                                  m[group==groups[i]],
                                  n[group==groups[i]],
                                  pt.method=pt.method,
                                  ci.method=ci.method,
                                  scale=scale,alpha=alpha,tol=tol)
        class(ans[[i]]) <- "pooledBin"
        #attributes(ans[[i]]) <- list(class="pooledBin",call=call)
        attr(ans[[i]],"call") <- call

      }
      names(ans) <- groups
      ans.lst <- ans
      ans <- data.frame(Group = attr(x, "group.names"),
                        PointEst = rep(0,  nGroups),
                        Lower = rep(0, nGroups),
                        Upper = rep(0, nGroups),
                        Scale = rep(1,  nGroups))
      if (!is.null(attr(ans.lst, "group.var")))
        names(ans)[1] <- attr(ans.lst, "group.var")
      for (i in 1:n) out[i, 2:5] <- ans.lst[[i]]$scale * c(ans.lst[[i]]$p,  ans.lst[[i]]$lcl, ans.lst[[i]]$ucl, 1)
      if (all(ans$Scale == 1)) ans$Scale <- NULL

      ans <- structure(ans, class = "pooledBinList", fullList = ans.lst,
                       x.var = vars$x, m.var = vars$m, n.var = vars$n,
                       group.names = groups, group.var = group.var,scale=scale,call=call)


    } else {
      ans <- pooledBin.fit(x,m,n,
                           pt.method=pt.method,
                           ci.method=ci.method,
                           scale=scale,alpha=alpha,tol=tol)
      #names(ans) <- "1"
      ans.lst <- ans
      ans <- data.frame(P=scale*ans.lst$p,
                        Lower = scale*ans.lst$lcl,
                        Upper = scale*ans.lst$ucl,
                        Scale = scale) # really to match Hmisc's binconf()
      if(scale == 1) ans$Scale <- NULL
      ans <- structure(ans,class="pooledBin",fullList = ans.lst,call=call,group.names = "", group.var = "", scale=scale)
      #ans <- structure(ans, class = "pooledBin", group.names = names(ans),
      #                x.var = vars$x, m.var = vars$m, n.var = vars$n, group.var = vars$group, scale=scale,call = call)

    }
    #attributes(ans,"x.var") <- vars$x
    #attributes(ans,"m.var") <- vars$m
    #attributes(ans,"n.var") <- vars$n
    #if(!is.null(vars$group)) attributes(ans,"group.var") <- vars$group

    ans

  }

"pooledBinParseFormula" <- function(formula,data=NULL){
  # this function reads a formula of the form
  # x ~ m(m.var) + n(n.var) | group
  # for a call to pooledBin and returns
  # a list with elements (x, m, n, group), which are character strings
  # with the respective names for those variables;
  # missing n() or group variables are set to NULL;
  # if a single variable is used it is assumed to be for m(), in which case
  # the m() construct may be excluded;
  # a grouping variable is optional;
  # errors for ambiguous calls are issued;

  # first determine if there's a grouping variable, indicated by "|" in the formula
  # note: grepl returns logicals for inclusion, whereas grep does not
  has.group <- any(grepl("\\|",formula))

  if(has.group){
    tmp.tms <- terms(formula,specials=c("m","n"),data=data)
    if(attr(tmp.tms,"response")==0) stop("must have a response variable representing the number of positives in the pool")
    group.var<- as.character(formula[[3]])[3]
    lhs <- as.character(formula)[2]
    rhs <- as.character(formula[[3]])[2]
    formula <- as.formula(paste0(lhs," ~ ",rhs))
  } else group.var <- NULL

  tms <- terms(formula, specials = c("m","n"),data=data)
  spec <- attr(tms,"specials")
  if(attr(tms,"response")==0) stop("must have a response variable representing the number of positives in the pool")
  if(length(attr(tms,"variables")) > 4) stop("at most pool size [m()] and number of pools [n()] accepted in RHS of formula...see help page for specification")
  x.var <- as.character(attr(tms,"variables")[[2]])
  m.ind <- spec$m
  if(!is.null(m.ind)){
    m.var <- as.character(attr(tms,"variables")[[m.ind+1]][[2]])
  } else m.var <- NULL
  n.ind <- spec$n
  # ifelse didn't want to play nicely with a NULL return
  if(!is.null(n.ind)){
    n.var <- as.character(attr(tms,"variables")[[n.ind+1]][[2]])
  } else n.var <- NULL

  if(is.null(m.var) & !is.null(n.var)){
    #stop("formula must contain an m() term if it contains an n() term\n")
    m.ind <- ifelse(n.ind==2,3,2)
    m.var <- as.character(attr(tms,"variables")[[m.ind+1]])
  }

  if(is.null(m.var) & is.null(n.var)) {
    if(length(attr(tms,"variables")) > 3) stop("if no m() or n() term specified, only one variable is permitted -- interpreted as the m() term")
    m.var <- as.character(attr(tms,"variables")[[3]])
  }

  list(x = x.var, m = m.var, n = n.var, group = group.var)
}


"pooledBin.fit" <- function(x,m,n=rep(1,length(x)),
           pt.method = c("firth","gart","bc-mle","mle","mir"),
           ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
           scale=1, alpha=0.05, tol=.Machine$double.eps^0.5,...)
  {
  # check for x > n
  if(any(x>n)) stop("must have x <= n for each element")
  # if all m==1, there is no pooling, so just use the regular, Wilson score interval
  if(all(m==1)){
    z <- qnorm(1-alpha/2)
    p <- sum(x)/sum(n)
    # the usual score interval, when all pools are size 1
    cl <- (p + z^2/2/sum(n) + c(-1, 1) * z * sqrt((p * (1 -  p) + z^2/4/sum(n))/sum(n)))/(1 + z^2/sum(n))

    return(list(p=p,lcl=cl[1],ucl=cl[2],pt.method="mle",ci.method="score",alpha=alpha,x=x,m=m,n=n,scale=scale))
  }
  # check for all x == n
  if(all(x == n)){
    #warning("estimation from pooled samples with all pools positive not available--estimates set to NA")
    return(list(p=NA,lcl=NA,ucl=NA,pt.method=pt.method,ci.method=ci.method,alpha=alpha,x=x,m=m,n=n,scale=scale))
  }

    pt.method <- match.arg(pt.method)
    ci.method <- match.arg(ci.method)
    if(ci.method=="mir" | pt.method=="mir"){
      ci.method <- "mir"
      pt.method <- "mir"
    }
    if(pt.method == "gart") pt.method <- "bc-mle" # backward compatability
    switch(pt.method,
           "firth" = {p <- pooledbinom.firth(x,m,n,tol)},
           "mle" = { p <- pooledbinom.mle(x,m,n,tol)},
           "bc-mle" ={ p <- pooledbinom.cmle(x,m,n,tol)},
           "mir" = {p <- pooledbinom.mir(x,m,n)}
    )
    if(p < 0 & pt.method=="bc-mle"){
      pt.method <- "mle"
      warning("Bias-correction results in negative point estimate; using MLE\n")
      p <- pooledbinom.mle(x,m,n,tol)
    }
    switch(ci.method,
           "skew-score" = { ci.p <- pooledbinom.cscore.ci(x,m,n,tol,alpha)[2:3]},
           "bc-skew-score" ={ ci.p <- pooledbinom.bcscore.ci(x,m,n,tol,alpha)[2:3]},
           "score" = {ci.p <- pooledbinom.score.ci(x,m,n,tol,alpha)[2:3]},
           "lrt" = {ci.p <- pooledbinom.lrt.ci(x,m,n,tol,alpha)[2:3]},
           "wald" = {ci.p <- pooledbinom.wald.ci(x,m,n,tol,alpha)[2:3]},
           "mir" = {ci.p <- pooledbinom.mir.ci(x,m,n,tol,alpha)[2:3]}
    )
    list(p=p,lcl=ci.p[1],ucl=ci.p[2],pt.method=pt.method,ci.method=ci.method,alpha=alpha,x=x,m=m,n=n,scale=scale)
  }


"print.pooledBin" <- function(x, ...){
    args <- list(...)
    #if(is.null(args$digits)) digits <- 4
    #else digits <- args$digits
    #p <- round(scale*x$p,digits)
    #lcl <- round(scale*x$lcl, digits)
    #ucl <- round(scale*x$ucl, digits)
    # p <- scale*x$p
    # lcl <- scale*x$lcl
    # ucl <- scale*x$ucl
    # mat <- matrix(c(p,lcl,ucl,scale),nrow=1) # really to match Hmisc's binconf()
    # dimnames(mat) <- list(c(""),c("P","Lower","Upper","Scale"))
    # if(scale == 1) mat <- mat[,-4]
    # print(mat,...)
    x
  }

"print.pooledBinList" <- function (x, ...)
{
  # n <- length(x)
  # out <- data.frame(Group = attr(x, "group.names"),
  #                   PointEst = rep(0,  n),
  #                   Lower = rep(0, n),
  #                   Upper = rep(0, n),
  #                   Scale = rep(1,  n))
  # if (!is.null(attr(x, "group.var")))
  #   names(out)[1] <- attr(x, "group.var")
  # for (i in 1:n) out[i, 2:5] <- x[[i]]$scale * c(x[[i]]$p,  x[[i]]$lcl, x[[i]]$ucl, 1)
  # if (all(out$Scale == 1))
  #   out$Scale <- NULL
  x
}

# "as.data.frame.pooledBin" <- function(x, row.names = NULL, optional = FALSE, ...){
#   args <- list(...)
#   #if(is.null(args$digits)) digits <- 4
#   #else digits <- args$digits
#   scale <- x$scale
#   #p <- round(scale*x$p,digits)
#   #lcl <- round(scale*x$lcl, digits)
#   #ucl <- round(scale*x$ucl, digits)
#   p <- scale*x$p
#   lcl <- scale*x$lcl
#   ucl <- scale*x$ucl
#   df <- data.frame(P=p,Lower=lcl,Upper=ucl,Scale=scale)
#   rownames(df) <- 1:nrow(df)
#   if(scale == 1) df$Scale <- NULL
#   df
# }

# "as.data.frame.pooledBinList" <- function (x, row.names = NULL, optional = FALSE, ...)
# {
#   n <- length(x)
#   out <- data.frame(Group = attr(x, "group.names"),
#                     PointEst = rep(0,  n),
#                     Lower = rep(0, n),
#                     Upper = rep(0, n),
#                     Scale = rep(1,  n))
#   if (!is.null(attr(x, "group.var")))
#     names(out)[1] <- attr(x, "group.var")
#   for (i in 1:n) out[i, 2:5] <- x[[i]]$scale * c(x[[i]]$p,  x[[i]]$lcl, x[[i]]$ucl, 1)
#   if (all(out$Scale == 1))
#     out$Scale <- NULL
#   out
# }

"summary.pooledBin" <-
  function(object, ...){
    x <- attr(object,"fullList")
    args <- list(...)
    scale <- x$scale
    #if(is.null(args$digits)) digits <- 4
    #else digits <- args$digits
    switch(x$pt.method,
           "firth" = x$PtEstName <- "Firth's Correction",
           "gart" = x$PtEstName <- "Gart's Correction",
           "bc-mle" = x$PtEstName <- "Gart's Correction",
           "mle" = x$PtEstName <- "Maximum Likelihood",
           "mir" = x$PtEstName <- "Minimum Infection Rate"
    )
    switch(x$ci.method,
           "skew-score" = x$CIEstName <- "Skew-Corrected Score (Gart)",
           "bc-skew-score" = x$CIEstName <- "Bias- & Skew-Corrected Score (Gart)",
           "score" = x$CIEstName <- "Score",
           "lrt" = x$CIEstName <- "Likelihood Ratio Test Inversion",
           "wald" = x$CIEstName <- "Wald",
           "mir" = x$CIEstName <- "Minimum Infection Rate"
    )
    structure(x,class="summary.pooledBin",call=attr(object,"call"),scale=x$scale)
  }

"print.summary.pooledBin" <-
  function(x, ...){
    args <- list(...)
    scale <- x$scale
    cat("Estimation of Binomial Proportion for Pooled Data\n\n")
    print.pooledBin(x, ...)
    cat("\n")
    cat(paste0("\nCall: ", deparse(attr(x,"call"),width.cutoff = 100),"\n\n"))
    cat(paste("Point estimator:",x$PtEstName,"\n"))
    cat(paste("CI method:",x$CIEstName,"\n\n"))
    cat(paste("Number of individuals:",sum(x$n * x$m),"\n"))
    cat(paste("Number of pools:",sum(x$n),"\n"))
    cat(paste("Number of positive pools:",sum(x$x),"\n"))
    invisible(x)
  }


"summary.pooledBinList" <- function(object, ...){
  grp.names <- as.character(attr(object,"group.names"))
  "sumf" <- function(x) c(x, N = sum(x$n * x$m), NumPools = sum(x$n), NumPosPools = sum(x$x),
                          PtEstName = x$pt.method,
                          CIEstName = x$ci.method,
                          Alpha = x$alpha) # c() just adds to the list x
  out <- lapply(attr(object,"fullList"), sumf)
  #attributes(out) <- list(class = "summary.pooledBinList", names = attr(object,"names"),group.var = attr(object,"group.var"))
  #attributes(out) <- list(class = "summary.pooledBinList",
  #                        names = grp.names,group.var = attr(object,"group.var"))
  structure(out, class = "summary.pooledBinList",
            names = grp.names,group.var = attr(object,"group.var"),
            call = attr(object,"call"))
}

"print.summary.pooledBinList" <- function(x, ...){
  n <- length(x)
  out <- data.frame(Group = names(x),
                    PointEst = rep(0,n),
                    Lower = rep(0,n),
                    Upper = rep(0,n),
                    N = rep(0,n),
                    NumPools = rep(0,n),
                    NumPosPools = rep(0,n),
                    Scale = rep(1,n))
  if(!is.null(attr(x,"group.var"))) names(out)[1] <- attr(x,"group.var")
  for(i in 1:n)
    out[i,2:8] <- c(x[[i]]$scale * c(x[[i]]$p, x[[i]]$lcl,x[[i]]$ucl), x[[i]]$N, x[[i]]$NumPools, x[[i]]$NumPosPools, x[[i]]$scale)

  #if(is.null(digits)) digits <- 4
  #else digits <- args$digits
  cat("\nEstimation of Binomial Proportion for Pooled Data\n")
  cat(paste0("\nCall: ", deparse(attr(x,"call"),width.cutoff = 100),"\n\n"))
  switch(x[[1]]$PtEstName,
         "firth"  = PtEstName <- "Firth's Correction",
         "gart"   = PtEstName <- "Gart's Correction",
         "bc-mle" = PtEstName <- "Gart's Correction",
         "mle"    = PtEstName <- "Maximum Likelihood",
         "mir"    = PtEstName <- "Minimum Infection Rate"
  )
  switch(x[[1]]$CIEstName,
         "skew-score"    = CIEstName <- "Skew-Corrected Score (Gart)",
         "bc-skew-score" = CIEstName <- "Bias- & Skew-Corrected Score (Gart)",
         "score"         = CIEstName <- "Score",
         "lrt"           = CIEstName <- "Likelihood Ratio Test Inversion",
         "wald"          = CIEstName <- "Wald",
         "mir"           = CIEstName <- "Minimum Infection Rate"
  )
  cat(paste("Point estimator        :",PtEstName,"\n"))
  cat(paste("CI method              :",CIEstName,"\n"))
  cat(paste("Confidence coefficient : ",100*(1-x[[1]]$alpha),"%\n\n",sep=""))

  print(out,...)
  invisible(x)
}




"plot.pooledBin" <-
  function(x,pch=16,refline=TRUE,printR2=TRUE,...){
    # reference: Chen & Swallow
    x.lst <- attr(x,"fullList")
    if(all(x.lst$n==1)) {
      xmn <- as.list(by(x.lst$x,x.lst$m,function(u) c(sum(u),length(u))))
      m <- as.numeric(names(xmn))
      xx <- sapply(xmn,function(u) u[1])
      n <- sapply(xmn,function(u) u[2])
    } else {
      xx <- x.lst$x
      m <- x.lst$m
      n <- x.lst$n
    }
    y <- log((xx+0.5)/(n+0.5))
    cc <- lm(y ~ m)
    plot(m,y,...)
    if(refline) abline(cc)
    if(printR2) cat(paste("R-squared for diagnostic line fit =",round(summary(cc)$r.squared,4),"\n"))
    invisible(x)
  }

"plot.pooledBinList" <-
  function(x,pch=16,refline=TRUE,printR2=TRUE,layout=NULL,...){
    x.lst <- attr(x,"fullList")
    n.groups <- length(x.lst)
    n.groups.root <- ceiling(sqrt(n.groups))
    if(is.null(layout)) layout <- c(n.groups.root, n.groups.root)
    par(mfrow = layout)
    for(i in 1:n.groups){
      # reference: Chen & Swallow
      if(all(x.lst[[i]]$n==1)) {
        xmn <- as.list(by(x.lst[[i]]$x, x.lst[[i]]$m, function(x) c(sum(x),length(x))))
        m <- as.numeric(names(xmn))
        xx <- sapply(xmn,function(x) x[1])
        n <- sapply(xmn,function(x) x[2])
      } else {
        xx <- x.lst[[i]]$x
        m <- x.lst[[i]]$m
        n <- x.lst[[i]]$n
      }
      y <- log((xx+0.5)/(n+0.5))
      cc <- lm(y ~ m)
      plot(m,y,xlab="Pool Size",ylab="log((x+0.5)/(n+0.5))",...)
      if(refline) abline(cc)
      #if(printR2) cat(paste("R-squared for diagnostic line fit =",round(summary(cc)$r.squared,4),"\n"))
      #if(printR2) title(expression(paste(names(x)[[i]],"\n(",plain(R)^2,"=",round(summary(cc)$r.squared,2),")")),cex=0.75)
      #if(printR2) title(bquote(atop(.(names(x)[[i]]),"(" ~ R^2 ~ "=" ~ .(round(summary(cc)$r.squared,2)) ~ ")")), cex = 0.75)
      if(printR2) title(bquote(.(names(x)[[i]]) ~ ":" ~ R^2 ~ "=" ~ .(round(summary(cc)$r.squared,2))), cex = 0.75)
      else title(names(x)[[i]],cex=0.75)
    }
    invisible(x)
  }
