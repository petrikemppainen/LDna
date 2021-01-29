#' Produces necessary data for LDna from a matrix of pairwise LD values.
#'
#' Takes a lower diagonal matrix of pairwise LD values and produces data files for subsequent LD network analyses. This version works in conjunction with \code{\link{extractBranches}}.
#'
#' @details In LDna clusters represent loci (vertices) connected by LD values (edges) above given thresholds. From a lower diagonal matrix of LD values a \code{a clustering algorithim} (see \code{\link{hclust}}) is used to produce a tree that describes cluster merging with decreasing LD threshold. In these trees nodes represent clusters/loci and node distance indicates LD threshold values at which these events occur. This version (in contrast to \code{\link{LDnaRaw}}) does not estimate lambda as it is intended to be used with \code{\link{extractBranches}} where |E|min (minimum number of edges) is the only parameter that defines LD-clusters (essentially determines how "thick" branches are).
#' \cr
#' \cr
#' The clustering method \code{method} can be specified. For 'strict' clusters use 'single' (e.g. for outlier or inversion detection), for chromosome binning, use 'average', 'complete' or 'median'.
#' R2 values can be rounded to limit (\code{digits}) number of clustering events in large data sets, but too much leads to problems when recursively collapsing (too large) polytomies.
#'
#' @param LDmat Lower diagonal matrix of pairwise LD values, \eqn{R^2} is strongly recommended. Must include row and column names.
#' @param digits Number of digits for rounding r2 values. 
#' @param method Specifies clustering medhod (see \code{\link{hclust}}), defaults to 'single'.
#' @param mc.cores Specifies number of cores for \code{\link{mclapply}}. \code{NULL} (default) uses \code{\link{lapply}}.
#' @param length_out Specifies the number of maximum number of possible node depths. Overrides \code{digits}, further details are given below.
#' @param fun Only used internally
#' @keywords LDnaRaw2
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, Christopher Knight \email{Chris.Knight@@manchester.ac.uk}
#' 
#' @return A list with three objects: \code{clusterfile}, \code{stats} and a tree-file \code{tree}. \code{Clusterfile} is a matrix with all unique clusters as columns and loci as rows, where \code{TRUE} indicates the presence of a locus in a specific cluster (else \code{FALSE} is given). \code{Stats} is a data frame that gives all edges for the \code{single linkage clustering tree} along with the following information for each edge: cluster (focal cluster), parent cluster (cluster efter merger), and the corresponding \code{nV} (number of loci), \code{nE}(number of edges connecting the loci), \code{min_LD} (minimum LD among any two loci in the cluster), \code{merge_at_below} (the weakest link in the cluster before merger at lower thresholds), \code{merge_at_above} (the weakest link in the cluster just before merger spliting into to or more clusters, at higher thresholds). See example below
#' @seealso \code{\link{extractBranches}}
#' @examples
#' # Simple upper diagonal LD matrix
#' LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.44, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.36, 0.24, NA, NA, NA, 0.48, 0.32, 0.2, 0.36, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.24, NA, NA, NA, NA, NA, NA, 0.44, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
#' # Produce raw data
#' ldna <- LDnaRaw2(LDmat)
#' ldna$clusterfile
#' ldna$stats
#' 
#' 
#' # Illustrating \code{merge_at_below} and \code{merge_at_above} 
#' clusters <- extractBranches(ldna, min.edges=0, merge.min=0.8)
#' clusters
#' cl_info <- ldna$stats[cluster==names(clusters)[2]] # focus on cluster 5_0.68
#' cl_info
#' abline(v=cl_info$merge_at_below, col="blue")
#' abline(v=cl_info$merge_at_above, col="green")
#' 
#' 
#' #' # using multiple cores for a larger data set
#' data(LDna)
#' ldna <- LDnaRaw2(r2.baimaii_subs, mc.cores=4)
#' @export

LDnaRaw2 <- function(LDmat, digits=2, method='single', mc.cores=NULL, length_out=NULL, fun=function(x){min(x, na.rm=TRUE)}){

  if(!is.null(length_out)){
    tmp <- seq(0,1,length.out = length_out+1)
    
    LDmat[lower.tri(LDmat)] <- unlist(mclapply(1:nrow(LDmat), function(x){
      lapply(x:nrow(LDmat), function(y){
        tmp[which.min(abs(tmp-LDmat[y,x]))]
      } )
    }, mc.cores = cores))
    digits <- 10
  }
  
  if(is.na(LDmat[2,1])) LDmat <- t(LDmat)
  LDmat[LDmat<0] <- 0
  LDmat <- round(LDmat, digits)
  tree <- ape::as.phylo(hclust(as.dist(1-LDmat), method=method))
  tree <- ape::di2multi(tree)
  tree$edge.length <- round(tree$edge.length, digits+1) 
  tree <- ape::di2multi(tree)
  
  Ntips <- length(tree$tip.label)
  
  if(!is.null(mc.cores)){
    out <- parallel::mclapply(Ntips:length(tree$edge.length)+1, function(x) stats.fun(tree, LDmat, x), mc.cores = mc.cores, mc.preschedule = TRUE)
  }else{
    out <- lapply(Ntips:length(tree$edge.length)+1, function(x) stats.fun(tree, LDmat, x))
  }
  
  
  nV <- sapply(out, function(x){x$nV})
  nE <- sapply(out, function(x){x$nE})
  min_LD <- sapply(out, function(x){x$min_LD})
  tot.d <- sapply(out, function(x){x$tot.d})
  clusterfile <- do.call('cbind', lapply(out, function(x){x$loci}))
  
  tot.d2 <- 1-2*tot.d
  temp <- paste(tree$Nnode+1-c(1:tree$Nnode), round(tot.d2, digits)[order(round(tot.d2, digits))], sep="_")
  c.names <- temp[order(order(round(tot.d2, digits)))]
  colnames(clusterfile) <- c.names
  rownames(clusterfile) <- tree$tip.label
  
  LDmat[upper.tri(LDmat)] <- t(LDmat)[upper.tri(LDmat)]
  
  names <- c(tree$tip.label, c.names)
  stats <- data.frame(nV, nE, min_LD)
  rownames(stats) <- c.names
  
  d <- cbind(cbind(names[tree$edge[,2]], names[tree$edge[,1]]), c(round(tree$edge.length*2, digits)))
  d <- rbind(c(d[1,2], "root", min(round(tot.d2, digits))), d)
  d <- data.frame(d , stats[match(d[,1], c.names), ])
  
  names(d) <- c("cluster", "parent_cluster", "distance", "nV", "nE", "min_LD")
  d$nV[is.na(d$nV)] <- 1
  d$nE[is.na(d$nE)] <- 0
  rownames(d) <- 1:nrow(d)
  d <- d[rev(order(d$nV)),]
  d$cluster <- as.character(d$cluster)
  d$parent_cluster <- as.character(d$parent_cluster)
 
  merge_at_below <- sapply(unique(d$cluster)[-1], function(cl){
    
    if(as.vector(na.omit(d$parent_cluster[d$cluster==cl]))!="root"){
      as.numeric(strsplit(d[d$cluster %chin%   as.vector(na.omit(d$parent_cluster[d$cluster==cl])),]$cluster, "_")[[1]][2])
    }else{
      0
    }
    
  })
  d <- data.table(d, merge_at_below=c(0,merge_at_below))

  sapply(unique(d$cluster)[-1], function(cl){
    tmp <- d[parent_cluster==cl,]
    if(any(tmp$merge_at_below+as.numeric(tmp$distance)!=1)){
      d[cluster==cl, merge_at_above := as.numeric(strsplit(tmp[tmp$merge_at_below+as.numeric(tmp$distance)!=1,]$cluster, "_")[[1]][2])  ]
    }else{
      if(all(tmp$merge_at_below+as.numeric(tmp$distance)==1) & nrow(tmp)>0){
        d[cluster==cl, merge_at_above := as.numeric(strsplit(cl, "_")[[1]][2])]
      }else{
        d[cluster==cl, merge_at_above := 1]
      }
    }
  })

  out <- list(clusterfile, d, tree)
  
  names(out) <- c("clusterfile", "stats", 'tree')
  return(out)
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
              loci=loci,
              min_LD=min(LDmat.temp, na.rm=TRUE))
}
