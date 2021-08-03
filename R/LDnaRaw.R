#' Produces necessary data for LDna from a matrix of pairwise LD values
#'
#' Takes a lower diagonal matrix of pairwise LD values and produces data files for subsequent LD network analyses.
#'
#' @details In LDna clusters represent loci (vertices) connected by LD values (edges) above given thresholds. From a lower diagonal matrix of LD values a \code{a clustering algorithim} (see \code{\link{hclust}}) is used to produce a tree that describes cluster merging with decreasing LD threshold. In these trees nodes represent clusters/loci and node distance indicates LD threshold values at which these events occur. To identify mergers of distinct clusters, the change in median LD before and after merger, \eqn{\lambda}, is calculated for all clusters in the tree. \eqn{\lambda} is defined as \eqn{(x_1-x_2)*n} where \eqn{x_1} is the median of all pairwise LD values within the focal cluster, \eqn{x_2} is the median of all pairwise values between loci within the focal cluster and between the loci in the focal cluster and any other loci it merges with, and \eqn{n} is the number of loci in the focal cluster. High values of \eqn{\lambda}, relative to the median of all values across the tree, indicate the merging of clusters with highly distinct LD signals. For each cluster the number of vertices, \code{nV}, the number of edges, \code{nE} are recorded. \code{LDnaRaw} includes  \code{\link{mclapply}} to enable the use of multiple cores. 
#' The clustering method \code{method} can now be specfied. For 'strict' clusters use 'single' (e.g. for outlier or inversion detection), for chromosome binning, use 'average', 'complete' or 'median'.
#' R2 values can be rounded to limit (\code{digits}) number of clustering events in large data sets, but too much leads to problems when recursively collapsing (too large) polytomies.
#'
#' @param LDmat Lower diagonal matrix of pairwise LD values, \eqn{R^2} is strongly recommended. Must include row and column names.
#' @param digits Number of digits for rounding r2 values. 
#' @param method Specifies clustering medhod (see \code{\link{hclust}}), defaults to 'single'.
#' @param mc.cores Specifies number of cores for \code{\link{mclapply}}. \code{NULL} (default) uses \code{\link{lapply}}.
#' @param fun Only used internally
#' @keywords LDnaRaw
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, Christopher Knight \email{Chris.Knight@@manchester.ac.uk}
#' @return A list with two objects: \code{clusterfile} and \code{stats}. \code{Clusterfile} is a matrix with all unique clusters as columns and loci as rows, where \code{TRUE} indicates the presense of a locus in a specific cluster (else \code{FALSE} is given). \code{Stats} is a data frame that gives all edges for the \code{single linkage clustering tree} along with the following information for each edge: cluster (focal cluster), parent cluster (cluster efter merger), and the corresponding \code{nV}, \code{nE}, and \eqn{\lambda} for each cluster (cluster). Each cluster is given a unique name which indicates the highest LD threshold at which it is present.
#' @seealso \code{\link{extractClusters}}, \code{\link{summaryLDna}} and \code{\link{plotLDnetwork}}
#' @examples
#' # Simple upper diagonal LD matrix
#' LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.44, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.36, 0.24, NA, NA, NA, 0.48, 0.32, 0.2, 0.36, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.24, NA, NA, NA, NA, NA, NA, 0.44, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
#' # Produce raw data
#' ldna <- LDnaRaw(LDmat)
#' # With multicore
#' ldna <- LDnaRaw(LDmat, mc.cores=4)
#' # Larger data set
#' data(LDna)
#' ldna <- LDnaRaw(r2.baimaii_subs)
#' # With multicore, in this case only slightly faster. 
#' ldna <- LDnaRaw(r2.baimaii_subs, mc.cores=4)
#' @export

LDnaRaw <- function(LDmat, digits=2, method='single', mc.cores=NULL, fun=function(x){min(x, na.rm=TRUE)}){
  if(is.na(LDmat[2,1])) LDmat <- t(LDmat)
  LDmat[LDmat<0] <- 0
  LDmat <- round(LDmat, digits)
  tree <- ape::as.phylo(hclust(as.dist(1-LDmat), method=method))
  tree <- ape::di2multi(tree)
  tree$edge.length <- round(tree$edge.length, digits = digits+1) 
  tree <- ape::di2multi(tree)
  
  Ntips <- length(tree$tip.label)
  
  if(!is.null(mc.cores)){
    out <- parallel::mclapply(Ntips:length(tree$edge.length)+1, function(x) stats.fun(tree, LDmat, x), mc.cores = mc.cores, mc.preschedule = TRUE)
  }else{
    out <- lapply(Ntips:length(tree$edge.length)+1, function(x) stats.fun(tree, LDmat, x))
  }
  
  
  nV <- sapply(out, function(x){x$nV})
  nE <- sapply(out, function(x){x$nE})
  tot.d <- sapply(out, function(x){x$tot.d})
  clusterfile <- do.call('cbind', lapply(out, function(x){x$loci}))
  
  tot.d2 <- 1-2*tot.d
  temp <- paste(tree$Nnode+1-c(1:tree$Nnode), round(tot.d2, digits)[order(round(tot.d2, digits))], sep="_")
  c.names <- temp[order(order(round(tot.d2, digits)))]
  colnames(clusterfile) <- c.names
  rownames(clusterfile) <- tree$tip.label
  
  LDmat[upper.tri(LDmat)] <- t(LDmat)[upper.tri(LDmat)]
  
  if(!is.null(mc.cores)){
    lambda <- unlist(parallel::mclapply(2:ncol(clusterfile), function(x) lambda.fun(tree, LDmat, clusterfile, Ntips, x), mc.cores = mc.cores, mc.preschedule = TRUE))
  }else{
    lambda <- sapply(2:ncol(clusterfile), function(x) lambda.fun(tree, LDmat, clusterfile, Ntips, x))}
  
  names <- c(tree$tip.label, c.names)
  stats <- data.frame(nV, nE, lambda=c(0,lambda))
  rownames(stats) <- c.names
  
  d <- cbind(cbind(names[tree$edge[,2]], names[tree$edge[,1]]), c(round(tree$edge.length*2, digits)))
  d <- rbind(c(d[1,2], "root", min(round(tot.d2, digits))), d)
  d <- data.frame(d , stats[match(d[,1], c.names), ])
  names(d) <- c("cluster", "parent_cluster", "distance", "nV", "nE", "lambda")
  d$nV[is.na(d$nV)] <- 1
  d$nE[is.na(d$nE)] <- 0
  d$lambda[is.na(d$lambda)] <- 0
  rownames(d) <- 1:nrow(d)
  d <- d[rev(order(d$nV)),]
  
  lambda_min <- data.table(do.call(rbind, lapply(2:ncol(clusterfile), function(x){
    p <- tree$edge[,1][tree$edge[,2] == Ntips+x]-Ntips
    temp <- LDmat[clusterfile[,x],clusterfile[,x]]
    temp <- temp[lower.tri(temp)]
    LD1 <- fun(temp[temp!=0])
    temp <- c(temp, LDmat[clusterfile[,p] + clusterfile[,x]==1,clusterfile[,x]])
    LD2 <- fun(temp[temp!=0])
    c(LD1,LD2)
  })))
  
  out <- list(clusterfile, d, tree, lambda_min)
  
  names(out) <- c("clusterfile", "stats", 'tree', 'lambda_min')
  return(out)
}

lambda.fun <- function(tree, LDmat, clusterfile, Ntips, x){
  p <- tree$edge[,1][tree$edge[,2] == Ntips+x]-Ntips
  temp <- LDmat[clusterfile[,x],clusterfile[,x]]
  temp <- temp[lower.tri(temp)]
  LD1 <- median(temp, na.rm=TRUE)
  temp <- c(temp, LDmat[clusterfile[,p] + clusterfile[,x]==1,clusterfile[,x]])
  LD2 <- median(temp, na.rm=TRUE)
  (LD1 - LD2)*length(which(clusterfile[,x]))/nrow(clusterfile)*1000 
}

stats.fun <- function(tree, LDmat, x){
  tree.temp <- ape::extract.clade(tree, x)
  loci <- colnames(LDmat) %in% tree.temp$tip.label
  LDmat.temp <- LDmat[loci, loci]
  
  tot.d <- tree$edge.length[which(tree$edge[,1]==x)[1]]
  temp <- x
  while(tree$edge[,2][tree$edge[,1]==temp][1]>x){
    temp <- tree$edge[,2][tree$edge[,1]==temp][1]
    tot.d <- tot.d+tree$edge.length[which(tree$edge[,1]==temp)[1]]
  }
  
  out <- list(nV=length(which(loci)), 
              nE=length(which(LDmat.temp >= 1-2*tot.d)),
              tot.d=tot.d,
              loci=loci)
}
