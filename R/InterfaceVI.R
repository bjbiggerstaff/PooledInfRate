
"vectorIndex" <- function(x,...){
  UseMethod("vectorIndex")
}
"VI" <- function(x,...){
  UseMethod("VI")
}


"VI.default" <- function(x,m,n=rep(1,length(x)), vector, trap.time=rep(1,length(x)), group,
                         n.use.traps = TRUE,
                         n.use.na = FALSE,
                         pt.method = c("firth","gart","bc-mle","mle","mir"),
                         ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
                         scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
  call <- match.call()
  call[[1]] <- as.name("VI")

  pt.method <- match.arg(pt.method)
  ci.method <- match.arg(ci.method)
  if(!missing(group)) group.var <- deparse(substitute(group))
  if(missing(vector)) stop("Must have vectors specified to use vectorIndex().\n")
  vectors.var <- deparse(substitute(vector))
  Vectors <- unique(vector)
  nVectors <- length(Vectors)
  if(missing(trap.time)){
    traptime.var <- "VectorTime"
  } else {
    traptime.var <- deparse(substitute(trap.time))
  }
  # note missing X but valid m and n, to be used for vector abundance
  x.with.na <- any(is.na(x)) # this is used as a switch below
  if(x.with.na){
    which.x.na <- is.na(x)
    x.na <- x[which.x.na]
    m.na <- m[which.x.na]
    n.na <- n[which.x.na]
    vector.na <- vector[which.x.na]
    trap.time.na <- trap.time[which.x.na]
    if(!missing(group)) group.na <- group[which.x.na]

    x <- x[!which.x.na]
    m <- m[!which.x.na]
    n <- n[!which.x.na]
    vector <- vector[!which.x.na]
    trap.time <- trap.time[!which.x.na]
    if(!missing(group)) group <- group[!which.x.na]
  }

  if(!missing(group)){
    groups <- unique(group)
    nGroups <- length(groups)
    ans <- vector(mode="list",length=nGroups)
    ans.p <- vector(mode="list",length=nGroups)
    names(ans.p) <- groups
    ans.nbar <- vector(mode="list",length=nGroups)
    names(ans.nbar) <- groups
    for(i in 1:nGroups){
      ans[[i]] <- 0
      ans.p[[i]] <- vector(mode="list",length=nVectors)
      names(ans.p[[i]]) <- Vectors
      ans.nbar[[i]] <- vector(mode="list",length=nVectors)
      names(ans.nbar[[i]]) <- Vectors
      for(j in 1:nVectors){
        sub <- (group==groups[i]) & (vector == Vectors[j])
        tmp.pb <- pooledBin.fit(x[sub], m[sub], n[sub],
                                pt.method=pt.method,
                                ci.method=ci.method,
                                scale=scale,alpha=alpha,tol=tol)
        ans.p[[i]][[j]] <- tmp.pb$p
        # include pools not tested (x=NA) in the estimate of abundance
        # but only if there are missing values in x
        if(x.with.na & n.use.na){
          if(n.use.traps){
            ans.nbar[[i]][[j]] <- (sum(m[sub]*n[sub]) +
                                     sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
              (sum(trap.time[sub]) + sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
          } else {
            ans.nbar[[i]][[j]] <-  (sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
              (sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
          }
        } else {
          ans.nbar[[i]][[j]] <- sum(m[sub]*n[sub])/  sum(trap.time[sub])
        }

        ans[[i]] <- ans[[i]] + ans.nbar[[i]][[j]] * ans.p[[i]][[j]]
      }
    }
    ans <- structure(ans,class = "VIList", group.names = groups, group.var = group.var,
                     vector = Vectors, vectors.var = vectors.var, traptime.var = traptime.var, #deparse(substitute(vectors)),
                     p = ans.p, n = ans.nbar,
                     n.use.na = n.use.na, n.use.traps = n.use.traps,
                     pt.method = pt.method, ci.method = ci.method,
                     call=call)
  } else {
    ans <- 0
    ans.p <- vector(mode="list",length=nVectors)
    names(ans.p) <- Vectors
    ans.nbar <- vector(mode="list",length=nVectors)
    names(ans.nbar) <- Vectors
    for(j in 1:nVectors){
      sub <- vector == Vectors[j]
      tmp.pb <- pooledBin.fit(x[sub], m[sub], n[sub],
                              pt.method=pt.method,
                              ci.method=ci.method,
                              scale=scale,alpha=alpha,tol=tol)
      ans.p[[j]] <- tmp.pb$p
      # include pools not tested (x=NA) in the estimate of abundance
      # but only if there are missing values in x
      if(x.with.na & n.use.na){
        if(n.use.traps){
          ans.nbar[[j]] <- (sum(m[sub]*n[sub]) + sum(m.na[vector.na==Vectors[j]] * n.na[vector.na==Vectors[j]])) /
                               (sum(trap.time[sub]) + sum(trap.time.na[vector.na==Vectors[j]]))
        } else {
          ans.nbar[[j]] <- (sum(m.na[vector.na==Vectors[j]] * n.na[vector.na==Vectors[j]])) /
            (sum(trap.time.na[vector.na==Vectors[j]]))
        }
      } else {
        ans.nbar[[j]] <- sum(m[sub]*n[sub]) / sum(trap.time[sub]) # number of each vector scaled by collection effort ("trap nights")
      }
      ans <- ans + ans.nbar[[j]] * ans.p[[j]]
    }
    ans <- structure(ans,class = "VI", vector = Vectors,
                     vectors.var = vectors.var, traptime.var = traptime.var, #deparse(substitute(vectors)),
                     p = ans.p, n = ans.nbar, n.use.na = n.use.na, n.use.traps = n.use.traps,
                     pt.method = pt.method, ci.method = ci.method,
                     call = call)
  }
  ans
}

"vectorIndex.default" <- function(x,m,n=rep(1,length(x)), vector, trap.time=rep(1,length(x)), group,
                                                  n.use.traps = TRUE,
                                                  n.use.na = FALSE,
                                                  pt.method = c("firth","gart","bc-mle","mle","mir"),
                                                  ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
                                                  scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
  call <- match.call()
  call[[1]] <- as.name("vectorIndex")
  pt.method <- match.arg(pt.method)
  ci.method <- match.arg(ci.method)
  if(!missing(group)) group.var <- deparse(substitute(group))
  if(missing(vector)) stop("Must have vectors specified to use vectorIndex().\n")
  vectors.var <- deparse(substitute(vector))
  Vectors <- unique(vector)
  nVectors <- length(Vectors)
  if(missing(trap.time)){
    traptime.var <- "VectorTime"
  } else {
    traptime.var <- deparse(substitute(trap.time))
  }
  # note missing X but valid m and n, to be used for vector abundance
  x.with.na <- any(is.na(x)) # this is used as a switch below
  if(x.with.na){
    which.x.na <- is.na(x)
    x.na <- x[which.x.na]
    m.na <- m[which.x.na]
    n.na <- n[which.x.na]
    vector.na <- vector[which.x.na]
    trap.time.na <- trap.time[which.x.na]
    if(!missing(group)) group.na <- group[which.x.na]

    x <- x[!which.x.na]
    m <- m[!which.x.na]
    n <- n[!which.x.na]
    vector <- vector[!which.x.na]
    trap.time <- trap.time[!which.x.na]
    if(!missing(group)) group <- group[!which.x.na]
  }

  if(!missing(group)){
    groups <- unique(group)
    nGroups <- length(groups)
    ans <- vector(mode="list",length=nGroups)
    ans.p <- vector(mode="list",length=nGroups)
    names(ans.p) <- groups
    ans.nbar <- vector(mode="list",length=nGroups)
    names(ans.nbar) <- groups
    for(i in 1:nGroups){
      ans[[i]] <- 0
      ans.p[[i]] <- vector(mode="list",length=nVectors)
      names(ans.p[[i]]) <- Vectors
      ans.nbar[[i]] <- vector(mode="list",length=nVectors)
      names(ans.nbar[[i]]) <- Vectors
      for(j in 1:nVectors){
        sub <- (group==groups[i]) & (vector == Vectors[j])
        tmp.pb <- pooledBin.fit(x[sub], m[sub], n[sub],
                                pt.method=pt.method,
                                ci.method=ci.method,
                                scale=scale,alpha=alpha,tol=tol)
        ans.p[[i]][[j]] <- tmp.pb$p
        # include pools not tested (x=NA) in the estimate of abundance
        # but only if there are missing values in x
        if(x.with.na & n.use.na){
          if(n.use.traps){
           ans.nbar[[i]][[j]] <- (sum(m[sub]*n[sub]) +
                                  sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
                                 (sum(trap.time[sub]) + sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
          } else {
           ans.nbar[[i]][[j]] <-  (sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
                                 (sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
          }
        } else {
          ans.nbar[[i]][[j]] <- sum(m[sub]*n[sub])/  sum(trap.time[sub])
        }

        ans[[i]] <- ans[[i]] + ans.nbar[[i]][[j]] * ans.p[[i]][[j]]
      }
    }
    ans <- structure(ans,class = "vectorIndexList", group.names = groups, group.var = group.var,
                     vector = Vectors, vectors.var = vectors.var, traptime.var = traptime.var, #deparse(substitute(vectors)),
                     p = ans.p, n = ans.nbar,
                     n.use.na = n.use.na, n.use.traps = n.use.traps,
                     pt.method=pt.method, ci.method=ci.method,
                     call=call)
  } else {
    ans <- 0
    ans.p <- vector(mode="list",length=nVectors)
    names(ans.p) <- Vectors
    ans.nbar <- vector(mode="list",length=nVectors)
    names(ans.nbar) <- Vectors
    for(j in 1:nVectors){
      sub <- vector == Vectors[j]
      tmp.pb <- pooledBin.fit(x[sub], m[sub], n[sub],
                                 pt.method=pt.method,
                                 ci.method=ci.method,
                                 scale=scale,alpha=alpha,tol=tol)
      ans.p[[j]] <- tmp.pb$p
      # include pools not tested (x=NA) in the estimate of abundance
      # but only if there are missing values in x
      if(x.with.na & n.use.na){
        if(n.use.traps){
          ans.nbar[[j]] <- (sum(m[sub]*n[sub]) +
                        sum(m.na[vector.na==Vectors[j]] * n.na[vector.na==Vectors[j]])) /
                       (sum(trap.time[sub]) + sum(trap.time.na[vector.na==Vectors[j]]))
        } else {
          ans.nbar[[j]] <- (sum(m.na[vector.na==Vectors[j]] * n.na[vector.na==Vectors[j]])) /
                       (sum(trap.time.na[vector.na==Vectors[j]]))
        }
      } else {
          ans.nbar[[j]] <- sum(m[sub]*n[sub]) / sum(trap.time[sub]) # number of each vector scaled by collection effort ("trap nights")
      }
      ans <- ans + ans.nbar[[j]] * ans.p[[j]]
    }
    ans <- structure(ans,class = "vectorIndex", vector = Vectors,
                     vectors.var = vectors.var, traptime.var = traptime.var, #deparse(substitute(vectors)),
                     p = ans.p, n = ans.nbar, n.use.na = n.use.na, n.use.traps = n.use.traps,
                     pt.method = pt.method, ci.method = ci.method,
                     call = call)
  }
  ans
}


"VI.formula" <- function(x, data,
                         n.use.traps = TRUE,
                         n.use.na = FALSE,
                         pt.method = c("firth","gart","bc-mle","mle","mir"),
                         ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
                         scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
  call <- match.call()
  call[[1]] <- as.name("VI")
  pt.method <- match.arg(pt.method)
  ci.method <- match.arg(ci.method)
  #print(call)
  #call[[1]] <- as.name(strsplit(as.character(call[[1]]),"\\.")[[1]][1])
  if(missing(data))
    data <- environment(x)


  vars <- VIParseFormula(x,data)
  if(any(sapply(vars,length)>1)) stop("only variable names permitted in formula; perhaps use the default call")

  if(is.na(vars$vector) | is.null(vars$vector) | vars$vector=="")
    stop("Must have vectors specified to use vectorIndex().\n")
  vector <- eval(parse(text=vars$vector),data)
  Vectors <- unique(vector)
  nVectors <- length(Vectors)
  x <- eval(parse(text=vars$x), data)
  m <- eval(parse(text=vars$m), data)
  if(!is.null(vars$n))
    n <- eval(parse(text=vars$n), data)
  else n <- rep(1,length(x)) # default n
  if(vars$traptime == ""){
    trap.time <- rep(1,length(x))
  } else{
    trap.time <- eval(parse(text=vars$traptime),data)
  }
  if(!is.null(vars$group)){
    group <- eval(parse(text=vars$group), data)
    groups <- unique(group)
    nGroups <- length(groups)
  }  else{
    group <- rep("SINGLE",length(x))
    groups <- unique(group)
    nGroups <- 1
  }

  # note missing X but valid m and n, to be used for vector abundance
  x.with.na <- any(is.na(x)) # this is used as a switch below
  if(x.with.na){
    which.x.na <- is.na(x)
    x.na <- x[which.x.na]
    m.na <- m[which.x.na]
    n.na <- n[which.x.na]
    vector.na <- vector[which.x.na]
    trap.time.na <- trap.time[which.x.na]
    if(!missing(group)) group.na <- group[which.x.na]

    x <- x[!which.x.na]
    m <- m[!which.x.na]
    n <- n[!which.x.na]
    vector <- vector[!which.x.na]
    trap.time <- trap.time[!which.x.na]
    if(!missing(group)) group <- group[!which.x.na]
  }

  if(!is.null(vars$group)){
    groups <- unique(group)
    nGroups <- length(groups)
    ans <- vector(mode="list",length=nGroups)
    ans.p <- vector(mode="list",length=nGroups)
    names(ans.p) <- groups
    ans.nbar <- vector(mode="list",length=nGroups)
    names(ans.nbar) <- groups
    for(i in 1:nGroups){
      ans[[i]] <- 0
      ans.p[[i]] <- vector(mode="list",length=nVectors)
      names(ans.p[[i]]) <- Vectors
      ans.nbar[[i]] <- vector(mode="list",length=nVectors)
      names(ans.nbar[[i]]) <- Vectors
      for(j in 1:nVectors){
        sub <- (group == groups[i]) & (vector == Vectors[j])
        tmp.pb <- pooledBin.fit(x[sub], m[sub], n[sub],
                                pt.method=pt.method,
                                ci.method=ci.method,
                                scale=scale,alpha=alpha,tol=tol)
        ans.p[[i]][[j]] <- tmp.pb$p
        ans.nbar[[i]][[j]] <- sum(m[sub] * n[sub]) / sum(trap.time[sub])# number of each vector
        # include pools not tested (x=NA) in the estimate of abundance
        # but only if there are missing values in x
        if(x.with.na & n.use.na){
          if(n.use.traps){
            ans.nbar[[i]][[j]] <- (sum(m[sub]*n[sub]) +
                                     sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
              (sum(trap.time[sub]) + sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
          } else {
            ans.nbar[[i]][[j]] <- (sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
              (sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
          }
        } else {
          ans.nbar[[i]][[j]] <- sum(m[sub]*n[sub])/  sum(trap.time[sub])
        }

        # number of each vector scaled by collection effort ("trap nights")
        ans[[i]] <- ans[[i]] + ans.nbar[[i]][[j]] * ans.p[[i]][[j]]
      }
    }
    ans <- structure(ans,class = "VIList", group.names = groups,
                     group.var = vars$group,
                     vector = Vectors, vectors.var = vars$vector, traptime.var=vars$traptime,
                     p = ans.p, n = ans.nbar,
                     n.use.na = n.use.na, n.use.traps = n.use.traps,
                     pt.method = pt.method, ci.method = ci.method,
                     call=call)
  } else {
    ans <- 0
    ans.p <- vector(mode="list",length=nVectors)
    names(ans.p) <- Vectors
    ans.nbar <- vector(mode="list",length=nVectors)
    for(j in 1:nVectors){
      sub <- vector == Vectors[j]

      tmp.pb <-  pooledBin.fit(x[sub], m[sub], n[sub],
                               pt.method=pt.method,
                               ci.method=ci.method,
                               scale=scale,alpha=alpha,tol=tol)
      ans.p[[j]] <-tmp.pb$p
      ans.nbar[[j]] <- sum(m[sub] * n[sub]) / sum(trap.time[sub])
      # include pools not tested (x=NA) in the estimate of abundance
      # but only if there are missing values in x
      if(x.with.na & n.use.na){
        if(n.use.traps){
          ans.nbar[[j]] <- (sum(m[sub]*n[sub]) +
                              sum(m.na[vector.na==Vectors[j]] * n.na[vector.na==Vectors[j]])) /
            (sum(trap.time[sub]) + sum(trap.time.na[vector.na==Vectors[j]]))
        } else {
          ans.nbar[[j]] <- (sum(m.na[vector.na==Vectors[j]] * n.na[vector.na==Vectors[j]])) /
            (sum(trap.time.na[vector.na==Vectors[j]]))
        }
      } else {
        ans.nbar[[j]] <- sum(m[sub]*n[sub]) / sum(trap.time[sub]) # number of each vector scaled by collection effort ("trap nights")
      }
      ans <- ans + ans.nbar[[j]] * ans.p[[j]]
    }
    ans <- structure(ans,class = "VI", group.names=NULL,group.var=NULL,
                     vector = Vectors, vectors.var = vars$vector, traptime.var=vars$traptime,
                     p = ans.p, n = ans.nbar, n.use.na = n.use.na, n.use.traps = n.use.traps,
                     pt.method = pt.method, ci.method = ci.method,
                     call = call)
  }
  ans
}


"vectorIndex.formula" <- function(x, data,
                                  n.use.traps = TRUE,
                                  n.use.na = FALSE,
                                                  pt.method = c("firth","gart","bc-mle","mle","mir"),
                                                  ci.method = c("skew-score","bc-skew-score","score","lrt","wald","mir"),
                                                  scale=1, alpha=0.05, tol=.Machine$double.eps^0.5, ...) {
  call <- match.call()
  call[[1]] <- as.name("vectorIndex")
  pt.method <- match.arg(pt.method)
  ci.method <- match.arg(ci.method)
  #print(call)
  #call[[1]] <- as.name(strsplit(as.character(call[[1]]),"\\.")[[1]][1])
  if(missing(data))
    data <- environment(x)


  vars <- VIParseFormula(x,data)
  if(any(sapply(vars,length)>1)) stop("only variable names permitted in formula; perhaps use the default call")

  if(is.na(vars$vector) | is.null(vars$vector) | vars$vector=="")
    stop("Must have vectors specified to use vectorIndex().\n")
  vector <- eval(parse(text=vars$vector),data)
  Vectors <- unique(vector)
  nVectors <- length(Vectors)
  x <- eval(parse(text=vars$x), data)
  m <- eval(parse(text=vars$m), data)
  if(!is.null(vars$n))
    n <- eval(parse(text=vars$n), data)
  else n <- rep(1,length(x)) # default n
  if(vars$traptime == ""){
    trap.time <- rep(1,length(x))
  } else{
    trap.time <- eval(parse(text=vars$traptime),data)
  }
  if(!is.null(vars$group)){
    group <- eval(parse(text=vars$group), data)
    groups <- unique(group)
    nGroups <- length(groups)
  }  else{
    group <- rep("SINGLE",length(x))
    groups <- unique(group)
    nGroups <- 1
  }

  # note missing X but valid m and n, to be used for vector abundance
  x.with.na <- any(is.na(x)) # this is used as a switch below
  if(x.with.na){
    which.x.na <- is.na(x)
    x.na <- x[which.x.na]
    m.na <- m[which.x.na]
    n.na <- n[which.x.na]
    vector.na <- vector[which.x.na]
    trap.time.na <- trap.time[which.x.na]
    if(!missing(group)) group.na <- group[which.x.na]

    x <- x[!which.x.na]
    m <- m[!which.x.na]
    n <- n[!which.x.na]
    vector <- vector[!which.x.na]
    trap.time <- trap.time[!which.x.na]
    if(!missing(group)) group <- group[!which.x.na]
  }

  if(!is.null(vars$group)){
    groups <- unique(group)
    nGroups <- length(groups)
    ans <- vector(mode="list",length=nGroups)
    ans.p <- vector(mode="list",length=nGroups)
    names(ans.p) <- groups
    ans.nbar <- vector(mode="list",length=nGroups)
    names(ans.nbar) <- groups
    for(i in 1:nGroups){
      ans[[i]] <- 0
      ans.p[[i]] <- vector(mode="list",length=nVectors)
      names(ans.p[[i]]) <- Vectors
      ans.nbar[[i]] <- vector(mode="list",length=nVectors)
      names(ans.nbar[[i]]) <- Vectors
      for(j in 1:nVectors){
        sub <- (group == groups[i]) & (vector == Vectors[j])
        tmp.pb <- pooledBin.fit(x[sub], m[sub], n[sub],
                               pt.method=pt.method,
                               ci.method=ci.method,
                               scale=scale,alpha=alpha,tol=tol)
        ans.p[[i]][[j]] <- tmp.pb$p
        ans.nbar[[i]][[j]] <- sum(m[sub] * n[sub]) / sum(trap.time[sub])# number of each vector
        # include pools not tested (x=NA) in the estimate of abundance
        # but only if there are missing values in x
        if(x.with.na & n.use.na){
          if(n.use.traps){
          ans.nbar[[i]][[j]] <- (sum(m[sub]*n[sub]) +
                                   sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
            (sum(trap.time[sub]) + sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
          } else {
            ans.nbar[[i]][[j]] <- (sum(m.na[group.na==groups[i] & vector.na==Vectors[j]] * n.na[group.na==groups[i] & vector.na==Vectors[j]])) /
                                   (sum(trap.time.na[group.na==groups[i] & vector.na==Vectors[j]]))
          }
        } else {
          ans.nbar[[i]][[j]] <- sum(m[sub]*n[sub])/  sum(trap.time[sub])
        }

        # number of each vector scaled by collection effort ("trap nights")
        ans[[i]] <- ans[[i]] + ans.nbar[[i]][[j]] * ans.p[[i]][[j]]
      }
    }
    ans <- structure(ans,class = "vectorIndexList", group.names = groups,
                     group.var = vars$group,
                     vector = Vectors, vectors.var = vars$vector, traptime.var=vars$traptime,
                     p = ans.p, n = ans.nbar,
                     n.use.na = n.use.na, n.use.traps = n.use.traps,
                     pt.method = pt.method, ci.method=ci.method,
                     call=call)
  } else {
    ans <- 0
    ans.p <- vector(mode="list",length=nVectors)
    names(ans.p) <- Vectors
    ans.nbar <- vector(mode="list",length=nVectors)
    for(j in 1:nVectors){
      sub <- vector == Vectors[j]

      tmp.pb <-  pooledBin.fit(x[sub], m[sub], n[sub],
                             pt.method=pt.method,
                             ci.method=ci.method,
                             scale=scale,alpha=alpha,tol=tol)
      ans.p[[j]] <-tmp.pb$p
      ans.nbar[[j]] <- sum(m[sub] * n[sub]) / sum(trap.time[sub])
      # include pools not tested (x=NA) in the estimate of abundance
      # but only if there are missing values in x
      if(x.with.na & n.use.na){
        if(n.use.traps){
        ans.nbar[[j]] <- (sum(m[sub]*n[sub]) +
                            sum(m.na[vector.na==Vectors[j]] * n.na[vector.na==Vectors[j]])) /
          (sum(trap.time[sub]) + sum(trap.time.na[vector.na==Vectors[j]]))
        } else {
        ans.nbar[[j]] <- (sum(m.na[vector.na==Vectors[j]] * n.na[vector.na==Vectors[j]])) /
          (sum(trap.time.na[vector.na==Vectors[j]]))
        }
      } else {
        ans.nbar[[j]] <- sum(m[sub]*n[sub]) / sum(trap.time[sub]) # number of each vector scaled by collection effort ("trap nights")
      }
      ans <- ans + ans.nbar[[j]] * ans.p[[j]]
    }
    ans <- structure(ans,class = "vectorIndex", group.names=NULL,group.var=NULL,
                     vector = Vectors, vectors.var = vars$vector, traptime.var=vars$traptime,
                     p = ans.p, n = ans.nbar, n.use.na = n.use.na, n.use.traps = n.use.traps,
                     pt.method = pt.method, ci.method = ci.method,
                     call = call)
  }
  ans
}




"VIParseFormula" <- function(formula,data=NULL){
  # this function reads a formula of the form
  # x ~ m(m.var) + n(n.var) | group / vector
  # or
  # x ~ m(m.var) + n(n.var) / vector
  # for a call to VI and returns
  # a list with elements (x, m, n, group, vector), which are character strings
  # with the respective names for those variables;
  # missing n() or group variables are set to NULL;
  # if a single variable is used it is assumed to be for m(), in which case
  # the m() construct may be excluded;
  # a grouping variable is optional;
  # errors for ambiguous calls are issued;

  # first determine if there's a grouping variable, indicated by "|" in the formula
  # note: grepl returns logicals for inclusion, whereas grep does not
  has.group <- any(grepl("\\|",formula))
  has.vectors <- any(grepl("/",formula))
  has.traptime <- any(grepl(":",formula))

  if(!has.vectors) stop("formula must specify the variable for the 'vectors' using / after the base formula")

  # pick off the bits before the / to get the vector & traptime variables
  vectors.var <- unlist(strsplit(as.character(formula),"/"))[[4]]
  if(has.traptime) {
    splt <- strsplit(vectors.var,":")
    vectors.var <- splt[[1]][[1]]
    traptime.var <- splt[[1]][[2]]
  } else {
    traptime.var <- ""
  }


  # strip off the bits after the / symbol
  # ...this leaves a formula that has the same structure as the pooledBin (or pIR) ones
  # ...so we can use that formula deparser to get the correct names
  formula <- as.formula(paste0(unlist(strsplit(as.character(formula), "/"))[-4][c(2,1,3)],collapse=" "))
  #
  lst <- pooledBinParseFormula(formula,data)
  # add the vectors variable names
  lst$vector <- vectors.var
  lst$traptime <- traptime.var
  lst
}

"print.VI" <- function(x, ...){
  print(c(VI = x), ...)
  invisible(x)
}

"print.vectorIndex" <- function(x, ...){
  args <- list(...)
  print(c(VI = x), ...)
  invisible(x)
}


"print.VIList" <- function(x, ...){
  args <- list(...)
  #cat("Vector Index\n")
  ans.dat <- data.frame(attr(x,"group.names"),VI = unlist(x))
  names(ans.dat)[1] <- attr(x,"group.var")
  print(ans.dat, ...)
  invisible(ans.dat)
}

"print.vectorIndexList" <- function(x, ...){
  args <- list(...)
  #cat("Vector Index\n")
  ans.dat <- data.frame(attr(x,"group.names"),VI = unlist(x))
  names(ans.dat)[1] <- attr(x,"group.var")
  print(ans.dat, ...)
  invisible(ans.dat)
}

"summary.vectorIndex" <- function(object, ...){
  x <- object
  args <- list(...)
  vi <- unlist(x)
  p <- attr(x,"p")
  n <- attr(x,"n")
  structure(list(vi=as.vector(vi),vector=attr(x,"vectors"),p=p,n=n,vi.obj=object),
              class="summary.vectorIndex")
}

"summary.VI" <- function(object, ...){
  x <- object
  args <- list(...)
  vi <- unlist(x)
  p <- attr(x,"p")
  n <- attr(x,"n")
  structure(list(vi=as.vector(vi),vector=attr(x,"vectors"),p=p,n=n,vi.obj=object),
              class="summary.VI")
}

"print.summary.vectorIndex" <- function(x, ...){
  args <- list(...)
  if(!is.null(args$digits)) digits <- args$digits
  else digits <- 4
  out.dat <- data.frame(attr(x$vi.obj,"vector"),
                        unlist(x$n),
                        unlist(x$p),
                        unlist(x$n) * unlist(x$p))
  names(out.dat) <- c(attr(x$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")

  cat(paste0("Vector Index = ",format(x$vi,digits=digits),"\n"))

  cat(paste("\nCall: ", deparse(attr(x$vi.obj,"call"),width.cutoff=100,nlines=1),"\n\n"))

  cat(paste0("Prevalence estimate method : ",attr(x$vi.obj,"pt.method"),"\n"))
  cat(paste0("Use non-tested individuals in abundance estimate : ",attr(x$vi.obj,"n.use.na"),"\n"))
  cat(paste0("Use tested pools' pool sizes in abundance estimate : ",attr(x$vi.obj,"n.use.traps"),"\n"))
  cat("\n")
  cat("Detail by vector:\n\n")
  print(out.dat, ...)
  invisible(x)
}


"print.summary.VI" <- function(x, ...){
  args <- list(...)
  if(!is.null(args$digits)) digits <- args$digits
  else digits <- 4
  out.dat <- data.frame(attr(x$vi.obj,"vector"),
                        unlist(x$n),
                        unlist(x$p),
                        unlist(x$n) * unlist(x$p))
  names(out.dat) <- c(attr(x$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")

  cat(paste0("Vector Index = ",format(x$vi,digits=digits),"\n"))

  cat(paste("\nCall: ", deparse(attr(x$vi.obj,"call"),width.cutoff=100,nlines=1),"\n\n"))

  cat(paste0("Prevalence estimate method : ",attr(x$vi.obj,"pt.method"),"\n"))
  cat(paste0("Use non-tested individuals in abundance estimate : ",attr(x$vi.obj,"n.use.na"),"\n"))
  cat(paste0("Use tested pools' pool sizes in abundance estimate : ",attr(x$vi.obj,"n.use.traps"),"\n"))
  cat("\n")
  cat("Detail by vector:\n\n")
  print(out.dat,...)
  invisible(x)
}

"summary.vectorIndexList" <- function(object, ...){
  x <- object
  args <- list(...)

  vi <- unlist(x)
  structure(list(vi=as.vector(vi),vector=attr(x,"vectors"),
                 p=attr(x,"p"),n=attr(x,"n"),nGroups = length(x), nVectors =length(attr(x,"vector")),vi.obj=object),
            class="summary.vectorIndexList")

}


"summary.VIList" <- function(object, ...){
  x <- object
  args <- list(...)

  vi <- unlist(x)
  structure(list(vi=as.vector(vi),vector=attr(x,"vectors"),
                 p=attr(x,"p"),n=attr(x,"n"),nGroups = length(x), nVectors =length(attr(x,"vector")),vi.obj=object),
            class="summary.VIList")

}


"print.summary.VIList" <- function(x, digits = 4, ...){
  args <- list(...)
  nGroups <- x$nGroups
  nVectors <- x$nVectors
  out.dat <- data.frame(rep(attr(x$vi.obj,"group.names"),each=nVectors),
                        rep(attr(x$vi.obj, "vector"),nGroups),
                        as.vector(unlist(x$n)),
                        as.vector(unlist(x$p)),
                        as.vector(unlist(x$n) * unlist(x$p)))
  names(out.dat) <- c(attr(x$vi.obj,"group.var"),attr(x$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")

  cat(paste("\nCall: ", deparse(attr(x$vi.obj,"call"),width.cutoff=100,nlines=1),"\n\n"))

  cat(paste0("Prevalence estimate method : ",attr(x$vi.obj,"pt.method"),"\n"))
  cat(paste0("Use non-tested individuals in abundance estimate : ",attr(x$vi.obj,"n.use.na"),"\n"))
  cat(paste0("Use tested pools' pool sizes in abundance estimate : ",attr(x$vi.obj,"n.use.traps"),"\n"))
  cat("\n")

  print(x$vi.obj) # print the VIs by group

  cat("\n")
  cat("Detail by group and vector:\n\n")
  print(out.dat,...)
  invisible(x)
}

"print.summary.vectorIndexList" <- function(x, digits = 4, ...){
  args <- list(...)
  nGroups <- x$nGroups
  nVectors <- x$nVectors
  out.dat <- data.frame(rep(attr(x$vi.obj,"group.names"),each=nVectors),
                        rep(attr(x$vi.obj, "vector"),nGroups),
                        as.vector(unlist(x$n)),
                        as.vector(unlist(x$p)),
                        as.vector(unlist(x$n) * unlist(x$p)))
  names(out.dat) <- c(attr(x$vi.obj,"group.var"),attr(x$vi.obj,"vectors.var"),"Avg N","P","(Avg N) * P")

  cat(paste("\nCall: ", deparse(attr(x$vi.obj,"call"),width.cutoff=100,nlines=1),"\n\n"))

  cat(paste0("Prevalence estimate method : ",attr(x$vi.obj,"pt.method"),"\n"))
  cat(paste0("Use non-tested individuals in abundance estimate : ",attr(x$vi.obj,"n.use.na"),"\n"))
  cat(paste0("Use tested pools' pool sizes in abundance estimate : ",attr(x$vi.obj,"n.use.traps"),"\n"))
  cat("\n")

  print(x$vi.obj,...) # print the VIs by group

  cat("\n")
  cat("Detail by group and vector:\n\n")
  print(out.dat,...)
  invisible(x)
}





