#' Modification of Emmax to be used for Genome Wide Association (GWA) analyses using Groups of SNPs connected by high LD.
#'
#' This function is not meant to be a part of LDna but is provided here to show an example of how output from \code{\link{LDnClustering}} can be used for GWA studies. Potentially there will be a separate R-package in the future that focuses on GWA. For now we refer to code provided with the reference, see below.
#' 
#' See reference
#'
#' @param Grp Index indicating which LD-cluster SNPs belong to. Produced by \code{\link{LDnClustering}}.
#' @param Y Phenotypic value
#' @param X Original SNP data used for by \code{\link{LDnClustering}} (matrix with SNPs as columns and individuals as rows)
#' @param K Kinship matrix, as produced e.g by \code{A.mat}, package \code{rrBLUP}
#' @param B Number of permutations. If \code{NULL}, permutation is not performed
#' @keywords Genome wide association analyses, GWAS, Linkage disequilibrium, network, clustering
#' @seealso \code{\link{LDnClustering}}
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, zitong.li \email{lizitong1985@@gmail.com}
#' @return See reference
#' @references 
#' Li, Z., Kemppainen, P., Rastas, P., Merila, J. Linkage disequilibrium clustering-based approach for association mapping with tightly linked genome-wide data. Accepted to Molecular Ecology Resources.
#' @examples
#' ## see \code{\link{LDnClustering}}
#' @export
emmax_group <- function(Y,X,Grp,K,B=NULL) {
  
  n<-length(Y)
  m<-ncol(X)
  p <- length(unique(Grp))
  
  
  stopifnot(ncol(K) == n)
  stopifnot(nrow(K) == n)
  stopifnot(nrow(X) == n)
  
  #INTERCEPT
  
  Xo<-rep(1,n)
  
  #K MATRIX NORMALISATION
  
  K_norm<-(n-1)/sum((diag(n)-matrix(1,n,n)/n)*K)*K
  
  #NULL MODEL
  
  null<-emma.REMLE(Y,as.matrix(Xo),K_norm)
  
  pseudoh<-null$vg/(null$vg+null$ve)
  
  #EMMAX
  
  M<-solve(chol(null$vg*K_norm+null$ve*diag(n)))
  Y_t<-crossprod(M,Y)
  Xo_t<-crossprod(M,Xo)
  
  RSS<-NULL
  
  df1<-NULL
  
  for (j in 1:p) {
    
    X_t<-crossprod(M,X[,Grp==j])
    
    df1[j] <- dim(X_t)[2]+1
    
    RSS[j]<- sum(lsfit(cbind(Xo_t,X_t),Y_t,intercept = FALSE)$residuals^2)
    rm(X_t)
    
  }
  
  
  RSSf<-unlist(RSS)
  RSS_H0<-rep(sum(lsfit(Xo_t,Y_t,intercept = FALSE)$residuals^2),p)
  df2<-n-df1-1
  R2<-1-1/(RSS_H0/RSSf)
  F<-(RSS_H0/RSSf-1)*df2/df1
  pval<-pf(F,df1,df2,lower.tail=FALSE)
  
  cat('EMMAX scan done! \n')
  
  
  ######  permutation
  
  if(!is.null(B)){
    
    sigma <- K_norm*null$vg+diag(rep(1,n))*null$ve
    #library(scalreg)
    Ynew <- mvrnorm(B, rep(0, ncol(sigma)), sigma)
    
    FR <- matrix(0,nrow=B,ncol=length(unique(Grp)))
    
    Fo <- sort(abs(F))
    
    Ord <- order(F)
    
    
    for (k in 1:B){
      
      y <- Ynew[k,]
      
      null<- emma.REMLE(y,as.matrix(Xo),K_norm)
      
      
      Y_t<-crossprod(M,y)
      
      for (j in 1:p) {
        
        X_t<-crossprod(M,X[,Grp==j])
        
        df1[j] <- dim(X_t)[2]+1
        
        RSS[j]<- sum(lsfit(cbind(Xo_t,X_t),Y_t,intercept = FALSE)$residuals^2)
        rm(X_t)
        
      }
      
      
      RSSf<-unlist(RSS)
      RSS_H0<-rep(sum(lsfit(Xo_t,Y_t,intercept = FALSE)$residuals^2),p)
      df2<-n-df1-1
      FR[k,] <-(RSS_H0/RSSf-1)*df2/df1
      
      
    }
    
    Qmat <- t(apply(FR,1,cummax))
    
    Padj <- apply(t(matrix(rep(Fo,B),length(unique(Grp)))) < Qmat, 2, mean)
    
    o <- order(Ord)
    
    return(list('F'=F,'pval'=pval,'pval.corr'=Padj[o],'Rsq'=R2))
    
  }else{
    return(pval)
  }
  
}

### not exported
emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10, esp=1e-10, eig.L = NULL, eig.R = NULL) {
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)
  
  #  stopifnot(nrow(K) == t)
  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)
  
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
    Etasq <- matrix(etas*etas,n-q,m)
    LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
    dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
      }
    }
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    Lambdas <- matrix(eig.R$values,t-q,m) + matrix(delta,t-q,m,byrow=TRUE)
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    dLL <- 0.5*delta*((n-q)*(colSums(Etasq/(Lambdas*Lambdas))+etas.2.sq/(delta*delta))/(colSums(Etasq/Lambdas)+etas.2.sq/delta)-(colSums(1/Lambdas)+(n-t)/delta))
    
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
      }
    }
  }  
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxva <- sum(etas*etas/(eig.R$values+maxdelta))/(n-q)    
  }
  else {
    maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/(n-q)
  }
  maxve <- maxva*maxdelta
  
  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
}

emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
  if ( complete == FALSE ) {
    vids <-  colSums(Z) > 0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  n <- nrow(Z)
  t <- ncol(Z)
  q <- ncol(X)
  
  SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
  eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE,EISPACK=TRUE)
  if ( is.complex(eig$values) ) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)    
  }
  qr.X <- qr.Q(qr(X))
  return(list(values=eig$values[1:(t-q)],
              vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                           complete=TRUE)[,c(1:(t-q),(t+1):n)]))   
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- lambda+delta
  return( 0.5*(nq*sum(etasq/(ldelta*ldelta))/sum(etasq/ldelta)-sum(1/ldelta)) )
}

emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(lambda+delta))))-sum(log(lambda+delta))) )
}


