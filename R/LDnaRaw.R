#' Produces necessary data for LDna from a matrix of pairwise LD values
#'
#' Takes a lower diagonal matrix of pairwise LD values and produces data files for subsequent LD network analyses.
#'
#' @details In LDna clusters represent loci (vertices) connected by LD values (edges) above given thresholds. From a lower diagonal matrix of LD values a \code{single linkage clustering algorithim} (see \code{\link{hclust}}) is used to produce a tree that describes cluster merging with decreasing LD threshold. In these trees nodes represent clusters/loci and node distance indicates LD threshold values at which these events occur. To identify mergers of distinct clusters, the change in median LD before and after merger, \eqn{\lambda}, is calculated for all clusters in the tree. \eqn{\lambda} is defined as \eqn{(x_1-x_2)*n} where \eqn{x_1} is the median of all pairwise LD values within the focal cluster, \eqn{x_2} is the median of all pairwise values between loci within the focal cluster and between the loci in the focal cluster and any other loci it merges with, and \eqn{n} is the number of loci in the focal cluster. High values of \eqn{\lambda}, relative to the median of all values across the tree, indicate the merging of clusters with highly distinct LD signals. For each cluster the number of vertices, \code{nV}, the number of edges, \code{nE} are recorded. \code{LDnaRaw} includes  \code{\link{mclapply}} to enable the use of multiple cores.  
#'
#' @param LDmat Lower diagonal matrix of pairwise LD values, \eqn{R^2} is strongly recommended.
#' @param mc.cores Specifies number of cores for \code{\link{mclapply}}. \code{NULL} (default) uses \code{\link{lapply}}.
#' @keywords LDnaRaw
#' @export
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
#' # With multicore, in this case only slightly faster
#' ldna <- LDnaRaw(r2.baimaii_subs, mc.cores=4)
# LDmat <- r2.baimaii_subs
LDnaRaw <- function(LDmat, mc.cores=NULL){
  LDmat1 <- round(LDmat, 2)
  temp2 <- hclust(as.dist(1-LDmat1), method="single")
  tree <- as.phylo(temp2)
  tree <- di2multi(tree)
  Ntips <- length(tree$tip.label)
  
  stats.fun <- function(x){
    tree.temp <- extract.clade(tree, x)
    loci <- colnames(LDmat1) %in% tree.temp$tip.label
    LDmat.temp <- LDmat1[loci, loci]
    
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
  
  if(!is.null(mc.cores)){out <- parallel::mclapply(Ntips:length(tree$edge.length)+1, stats.fun, mc.cores = mc.cores, mc.preschedule = TRUE)
  }else{out <- lapply(Ntips:length(tree$edge.length)+1, stats.fun)}
  
  
  nV <- sapply(out, function(x){x$nV})
  nE <- sapply(out, function(x){x$nE})
  tot.d <- sapply(out, function(x){x$tot.d})
  clusterfile <- do.call('cbind', lapply(out, function(x){x$loci}))
  
  tot.d2 <- 1-2*tot.d
  temp <- paste(tree$Nnode+1-c(1:tree$Nnode), round(tot.d2, 3)[order(round(tot.d2, 3))], sep="_")
  c.names <- temp[order(order(round(tot.d2, 3)))]
  colnames(clusterfile) <- c.names
  rownames(clusterfile) <- tree$tip.label
  
  LDmat3 <- LDmat
  LDmat3[upper.tri(LDmat3)] <- t(LDmat3)[upper.tri(LDmat3)]
  
  lambda.fun <- function(x){
    p <- tree$edge[,1][tree$edge[,2] == Ntips+x]-Ntips
    temp <- LDmat3[clusterfile[,x],clusterfile[,x]]
    temp <- temp[lower.tri(temp)]
    LD1 <- median(temp, na.rm=TRUE)
    temp <- c(temp, LDmat3[clusterfile[,p] + clusterfile[,x]==1,clusterfile[,x]])
    LD2 <- median(temp, na.rm=TRUE)
    list(lambda=(LD1 - LD2)*length(which(clusterfile[,x])))
  }
  
  if(!is.null(mc.cores)){out <- parallel::mclapply(2:ncol(clusterfile), lambda.fun, mc.cores = mc.cores, mc.preschedule = TRUE)
  }else{out <- lapply(2:ncol(clusterfile), lambda.fun)}
  
  lambda <- sapply(out, function(x){x$lambda})
  names <- c(tree$tip.label, c.names)
  stats <- data.frame(nV, nE, lambda=c(0,lambda))
  rownames(stats) <- c.names
  
  d <- tree$edge
  d[,2] <- names[tree$edge[,1]]
  d[,1] <- names[tree$edge[,2]]
  d <- cbind(d, c(round(tree$edge.length*2, 3)))
  d <- rbind(c(d[1,2], "root", min(round(tot.d2, 3))), d)
  d <- data.frame(d , stats[match(d[,1], c.names), ])
  names(d) <- c("cluster", "parent_cluster", "edge.length", "nV", "nE", "lambda")
  d$nV[is.na(d$nV)] <- 1
  d$nE[is.na(d$nE)] <- 0
  d$lambda[is.na(d$lambda)] <- 0
  d <- d[rev(order(d$nV)),]
  rownames(d) <- 1:nrow(d)
  out <- list(clusterfile, d)
  names(out) <- c("clusterfile", "stats")
  return(out)
}
