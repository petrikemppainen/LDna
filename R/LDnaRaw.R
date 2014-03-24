#' Produces necessary data for LDna from a matrix of pairwise LD values
#'
#' Takes a matrix of pairwise LD values and produces data files for subsequent LD network analyses.
#' 
#' In LDna, trees are used to visualise the process of cluster merging with decreasing LD threshold, where branches represent clusters/loci and branch points represent the merger of one or several clusters/loci. Distance indicates LD threshold values at which these events occur. To identify mergers of distinct clusters, the change in median LD between all loci in a cluster, before and after merger, lambda, is recorded for all clusters in the tree.  With a lower diagonal matrix of LD values, the LD values are transformed into dissimilarities (1-LD) and a \code{single linkage clustering} algorithim (see \code{\link{hclust}}) is used to produce a tree that describes the order of cluster merging with decreasing LD threshold. For each cluster the number of vertices, \code{nV}, the number of edges, \code{nE} are recorded.
#'
#' @param LDmat Lower diagonal matrix of pairwise LD values
#' @keywords LDnaRaw
#' @export
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, Christopher Knight \email{Chris.Knight@@manchester.ac.uk}
#' @return A list with two objects: \code{clusterfile} and \code{stats}. The first, \code{clusterfile} is a matrix with all unique clusters (defined by a specific set of loci) as columns and each locus as rows. In this matrix \code{TRUE} indicates the presense of a locus in a specific cluster, else \code{FALSE} is given. The second, \code{stats} is a data frame that gives all the edges for the \code{single linkage clustering tree} along with the following information for each edge: probe, parent probe, and the corresponding \code{nV}, \code{nE}, and the lambda for each probe.
#' @examples
#' # Simple upper diagonal LD matrix
#' LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.44, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.36, 0.24, NA, NA, NA, 0.48, 0.32, 0.2, 0.36, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.24, NA, NA, NA, NA, NA, NA, 0.44, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
#' # Calculate raw data
#' ldna <- LDnaRaw(LDmat)
#' # For larger data set
#' data(LDna)
#' ldna <- LDnaRaw(r2.baimaii_subs)

LDnaRaw <- function(LDmat){
  LDmat1 <- round(LDmat, 2)
  LDmat2 <- 1-LDmat1
  temp2 <- hclust(as.dist(LDmat2), method="single")
  phylo <- as.phylo(temp2)
  phylo <- di2multi(phylo)
  Ntips <- length(phylo$tip.label)
  
  nV <- rep(1, Ntips)
  nE <- rep(0, Ntips)
  tot.d <- NULL
  clusterfile <- matrix(FALSE, Ntips,  phylo$Nnode)
  for(i in Ntips:length(phylo$edge.length)+1){
    phylo.temp <- extract.clade(phylo, i)
    nV <- c(nV, length(phylo.temp$tip.label))
    LDmat.temp <- LDmat1[colnames(LDmat1) %in% phylo.temp$tip.label, rownames(LDmat1) %in% phylo.temp$tip.label]
    
    x <- phylo$edge.length[which(phylo$edge[,1]==i)[1]]
    temp <- i
    while(phylo$edge[,2][phylo$edge[,1]==temp][1]>i){
      temp <- phylo$edge[,2][phylo$edge[,1]==temp][1]
      x <- x+phylo$edge.length[which(phylo$edge[,1]==temp)[1]]
    }
    tot.d <- c(tot.d, x)
    nE <- c(nE, length(which(LDmat.temp >= 1-2*x)))
    
    temp <- i
    out <- i
    while(any(phylo$edge[,1] %in% temp)){
      temp <- phylo$edge[,2][phylo$edge[,1] %in% temp]
      out <- c(out, temp)
    }
    clusterfile[,i-Ntips][out[out<=Ntips]] <- TRUE
    
  }
  
  temp <- paste(phylo$Nnode+1-c(1:phylo$Nnode), round(1-2*tot.d, 3)[order(round(1-2*tot.d, 3))], sep="_")
  c.names <- temp[order(order(round(1-2*tot.d, 3)))]
  colnames(clusterfile) <- c.names
  rownames(clusterfile) <- phylo$tip.label

  Lambda <- rep(0, Ntips+phylo$Nnode)
  LDmat3 <- LDmat
  LDmat3[upper.tri(LDmat3)] <- t(LDmat3)[upper.tri(LDmat3)]
  for(i in 2:ncol(clusterfile)){
    #p <- phylo$edge[,1][phylo$edge[,2] == Ntips+i]-Ntips
    #LD1 <- median(LDmat3[clusterfile[,i],clusterfile[,i]], na.rm=T)
    #LD2 <- median(LDmat3[clusterfile[,p],clusterfile[,p]], na.rm=T)
    #Lambda[i+Ntips] <- (LD1 - LD2)*length(which(clusterfile[,i]))
    p <- phylo$edge[,1][phylo$edge[,2] == Ntips+i]-Ntips
    temp <- LDmat3[clusterfile[,i],clusterfile[,i]]
    temp <- temp[lower.tri(temp)]
    LD1 <- median(temp, na.rm=T)
    temp <- c(temp, LDmat3[clusterfile[,p] + clusterfile[,i]==1,clusterfile[,i]])
    LD2 <- median(temp, na.rm=T)
    Lambda[i+Ntips] <- (LD1 - LD2)*length(which(clusterfile[,i]))
  }
  
  stats <- data.frame(nV, nE, Lambda)
  rownames(stats) <- c(phylo$tip.label, c.names)
  
  d <- phylo$edge
  d[,2] <- rownames(stats)[phylo$edge[,1]]
  d[,1] <- rownames(stats)[phylo$edge[,2]]
  d <- cbind(d, c(round(phylo$edge.length*2, 3)))
  d <- rbind(c(d[1,2], "root", min(round(1-2*tot.d, 3))), d)
  d <- data.frame(d , stats[match(d[,1], rownames(stats)), ])
  names(d) <- c("probe", "parent_probe", "edge.length", "nV", "nE", "Lambda")
  d <- d[rev(order(d$nV)),]
  out <- list(clusterfile, d)
  names(out) <- c("clusterfile", "stats")
  return(out)
}

