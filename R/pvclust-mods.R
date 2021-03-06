####################################
# pvclust modifications (internal) #
####################################

text.pvclust <- function(x, col=c(2,3,8), print.num=TRUE,  float=0.01, cex=NULL, font=NULL, ...){
  axes <- hc2axes(x$hclust)
  usr  <- par()$usr; wid <- usr[4] - usr[3]
  au <- as.character(round(1-x$edges[,"au"],3))
  bp <- as.character(round(1-x$edges[,"bp"],3))
  rn <- as.character(row.names(x$edges))
  au[length(au)] <- "Approximately Unbiased (best choice)"
  bp[length(bp)] <- "Bootstrap Probability"
  rn[length(rn)] <- "edge #"
  a <- text(x=axes[,1], y=axes[,2] + float * wid, au,
            col=col[1], pos=2, offset=.3, cex=cex, font=font)
  a <- text(x=axes[,1], y=axes[,2] + float * wid, bp,col=col[2], pos=4, offset=.3, cex=cex, font=font)
  if(print.num)
    a <- text(x=axes[,1], y=axes[,2], "",col=col[3], pos=1, offset=.3, cex=cex, font=font)
}

hc2axes <- function(x){
  A <- x$merge # (n,n-1) matrix
  n <- nrow(A) + 1
  x.axis <- c()
  y.axis <- x$height
  
  x.tmp  <- rep(0,2)
  zz     <- match(1:length(x$order),x$order)
  
  for(i in 1:(n-1)) {
    ai <- A[i,1]
    if(ai < 0) x.tmp[1] <- zz[-ai] else x.tmp[1] <- x.axis[ai]
    
    ai <- A[i,2]
    if(ai < 0) x.tmp[2] <- zz[-ai] else x.tmp[2] <- x.axis[ai]
    
    x.axis[i] <- mean(x.tmp)
  }
  
  return(data.frame(x.axis=x.axis,y.axis=y.axis))
}

hc2split <- function(x){
  A <- x$merge # (n-1,n) matrix
  n <- nrow(A) + 1
  B <- list()
  
  for(i in 1:(n-1)){
    ai <- A[i,1]
    
    if(ai < 0) B[[i]] <- -ai
    else B[[i]] <- B[[ai]]        
    
    ai <- A[i,2]
    
    if(ai < 0) B[[i]] <- sort(c(B[[i]],-ai))
    else B[[i]] <- sort(c(B[[i]],B[[ai]]))
  }
  
  CC <- matrix(rep(0,n*(n-1)),nrow=(n-1),ncol=n)
  
  for(i in 1:(n-1)){
    bi <- B[[i]]
    m <- length(bi)
    for(j in 1:m) CC[i,bi[j]] <- 1
  }
  
  split <- list(pattern=apply(CC,1,paste,collapse=""), member=B)
  
  return(split)
}

pvclust.merge <- function(data, object.hclust, mboot){
  
  pattern <- hc2split(object.hclust)$pattern
  
  r     <- unlist(lapply(mboot,"[[","r"))
  nboot <- unlist(lapply(mboot,"[[","nboot"))
  store <- lapply(mboot,"[[", "store")
  
  rl <- length(mboot)
  ne <- length(pattern)
  
  edges.bp <- edges.cnt <- data.frame(matrix(rep(0,ne*rl),nrow=ne,ncol=rl))
  row.names(edges.bp) <- pattern
  names(edges.cnt) <- paste("r", 1:rl, sep="")
  
  for(j in 1:rl) {
    edges.cnt[,j] <- as.vector(mboot[[j]]$edges.cnt) 
    edges.bp[,j]  <- edges.cnt[,j] / nboot[j]
  }
  
  ms.fitted <- lapply(as.list(1:ne),
                      function(x, edges.bp, r, nboot){
                        msfit(as.vector(t(edges.bp[x,])), r, nboot)},
                      edges.bp, r, nboot)
  class(ms.fitted) <- "mslist"
  
  p    <- lapply(ms.fitted,"[[","p")
  se   <- lapply(ms.fitted,"[[","se")
  coef <- lapply(ms.fitted,"[[","coef")
  
  au    <- unlist(lapply(p,"[[","au"))
  bp    <- unlist(lapply(p,"[[","bp"))
  se.au <- unlist(lapply(se,"[[","au"))
  se.bp <- unlist(lapply(se,"[[","bp"))
  v     <- unlist(lapply(coef,"[[","v"))
  cc    <- unlist(lapply(coef,"[[","c"))
  pchi  <- unlist(lapply(ms.fitted,"[[","pchi"))
  
  edges.pv <- data.frame(au=au, bp=bp, se.au=se.au, se.bp=se.bp,
                         v=v, c=cc, pchi=pchi)
  
  row.names(edges.pv) <- row.names(edges.cnt) <- 1:ne
  
  result <- list(hclust=object.hclust, edges=edges.pv, count=edges.cnt,
                 msfit=ms.fitted, nboot=nboot, r=r, store=store)
  
  class(result) <- "pvclust"
  return(result)
}

dist.pvclust <- function(x, method="euclidean", use.cor="pairwise.complete.obs"){
  if(!is.na(pmatch(method,"correlation"))){
    res <- as.dist(1 - cor(x, method="pearson", use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(cor(x,method="pearson",use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  else if(!is.na(pmatch(method,"uncentered"))){
    if(sum(is.na(x)) > 0){
      x <- na.omit(x)
      warning("Rows including NAs were omitted")
    }
    x  <- as.matrix(x)
    P  <- crossprod(x)
    qq <- matrix(diag(P),ncol=ncol(P))
    Q  <- sqrt(crossprod(qq))
    res <- as.dist(1 - P/Q)
    attr(res,"method") <- "uncentered"
    return(res)
  }
  else
    dist(t(x),method)
}

corw <- function(x,w,use=c("all.obs","complete.obs","pairwise.complete.obs")){
  if(is.data.frame(x)) x <- as.matrix(x)
  x <- x[w>0,,drop=F]
  w <- w[w>0]
  
  n <- nrow(x) # sample size
  m <- ncol(x) # number of variables
  if(missing(w)) w <- rep(1,n)
  r <- matrix(0,m,m,dimnames=list(colnames(x),colnames(x)))
  diag(r) <- 1
  use <- match.arg(use)
  
  pairu <- F
  if(use=="all.obs") {
    u <- rep(T,n)
  } else if(use=="complete.obs") {
    u <- apply(x,1,function(y) !any(is.na(y)))
  } else if(use=="pairwise.complete.obs") {
    pairu <- T
    ux <- is.finite(x)
  } else stop("unknown use")
  
  for(i in 1+seq(length=m-1)) {
    for(j in seq(length=i-1)) {
      if(pairu) u <- ux[,i] & ux[,j]
      wu <- w[u]; xi <- x[u,i]; xj <- x[u,j]
      ws <- sum(wu)
      if(ws > 1e-8) {
        xi <- xi - sum(wu*xi)/ws
        xj <- xj - sum(wu*xj)/ws
        vxi <- sum(wu*xi*xi)/ws
        vxj <- sum(wu*xj*xj)/ws
        if(min(vxi,vxj) > 1e-8)  {
          vxij <- sum(wu*xi*xj)/ws
          rij <- vxij/sqrt(vxi*vxj)
        } else {
          rij <- 0
        }
      } else {
        rij <- 0
      }
      r[i,j] <- r[j,i] <- rij
    }
  }
  r
}

### calculate distance by weight
distw.pvclust <- function(x,w,method="correlation", use.cor="pairwise.complete.obs"){
  if(!is.na(pmatch(method,"correlation"))){
    res <- as.dist(1 - corw(x,w, use=use.cor))
    attr(res,"method") <- "correlation"
    return(res)
  }
  else if(!is.na(pmatch(method,"abscor"))){
    res <- as.dist(1 - abs(corw(x,w, use=use.cor)))
    attr(res,"method") <- "abscor"
    return(res)
  }
  stop("wrong method")
}


pvclust_new <- function(data, method.hclust="average",use.cor="pairwise.complete.obs",
                        nboot=1000, r=seq(.5,1.4,by=.1), store=FALSE, weight=FALSE){
  # data: (n,p) matrix, n-samples, p-variables
  n <- nrow(data); p <- ncol(data)
  
  # hclust for original data
  METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
               "median", "centroid")
  method.hclust <- METHODS[pmatch(method.hclust, METHODS)]
  distance <- mydist(data)
  data.hclust <- hclust(distance, method=method.hclust)
  
  # multiscale bootstrap
  size <- floor(n*r)
  rl <- length(size)
  
  if(rl == 1) {
    if(r != 1.0)
      warning("Relative sample size r is set to 1.0. AU p-values are not calculated\n")
    
    r <- list(1.0)
  }
  else
    r <- as.list(size/n)
  
  mboot <- lapply(r, boot.hclust_new, data=data, object.hclust=data.hclust, nboot=nboot,
                  use.cor=use.cor,
                  method.hclust=method.hclust, store=store, weight=weight)
  
  result <- pvclust.merge(data=data, object.hclust=data.hclust, mboot=mboot)
  
  return(result)
}

boot.hclust_new <- function(r, data, object.hclust,  use.cor,
                            method.hclust, nboot, store, weight=F){ 
  n     <- nrow(data)
  size  <- round(n*r, digits=0)
  if(size == 0)
    stop("invalid scale parameter(r)")
  r <- size/n
  
  pattern   <- hc2split(object.hclust)$pattern
  edges.cnt <- table(factor(pattern)) - table(factor(pattern))
  st <- list()
  
  # bootstrap start
  rp <- as.character(round(r,digits=2)); if(r == 1) rp <- paste(rp,".0",sep="")
  cat(paste("Bootstrap (r = ", rp, ")... ", sep=""))
  w0 <- rep(1,n) # equal weight
  na.flag <- 0
  
  for(i in 1:nboot){
    smpl <- sample(1:n, size, replace=TRUE)
    suppressWarnings(distance  <- mydist(data[smpl,]))
    
    if(all(is.finite(distance))) { # check if distance is valid
      x.hclust  <- hclust(distance,method=method.hclust)
      pattern.i <- hc2split(x.hclust)$pattern # split
      edges.cnt <- edges.cnt + table(factor(pattern.i,  levels=pattern))
    } else {
      x.hclust <- NULL
      na.flag <- 1
    }
    
    if(store)
      st[[i]] <- x.hclust
  }
  cat("Done.\n")
  # bootstrap done
  
  if(na.flag == 1)
    warning(paste("inappropriate distance matrices are omitted in computation: r = ", r), call.=FALSE)
  
  boot <- list(edges.cnt=edges.cnt, use.cor=use.cor,
               method.hclust=method.hclust, nboot=nboot, size=size, r=r, store=st)
  class(boot) <- "boot.hclust"
  
  return(boot)
}
