#' Summerises LD clusters
#'
#' Creates a data frame with summary information of LD clusters extracted by function \code{\link{extractClusters}}
#'
#' @param ldna Output from \code{\link{LDnaRaw}}
#' @param clusters Output from \code{\link{extractClusters}}
#' @param LDmat a matrix of pairwise LD values
#' @keywords summaryLDna
#' @seealso \code{\link{extractClusters}}, \code{\link{LDnaRaw}} and \code{\link{plotLDnetwork}}
#' @return Returns a data frame with each row corresponding to a cluster, in decreasing order with respect to highest LD threshold value at which they are present. Column \emph{Type} specifies if a cluster is a \emph{"single outlier cluster", SOC} or a \emph{"compound oulier cluster", COC}. Column \emph{"Merge.at"} specfies the LD threshold for cluster merger. "Median.LD" gives the average LD of all pairwise values between loci in a cluster, and "MAD.LD" gives their unscaled median absolute deviation.
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @examples
#' # Simple upper diagonal LD matrix
#' LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.44, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.36, 0.24, NA, NA, NA, 0.48, 0.32, 0.2, 0.36, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.24, NA, NA, NA, NA, NA, NA, 0.44, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
#' # Produce neccesary data for LDna
#' ldna <- LDnaRaw(LDmat)
#' # Extract clusters and plot graphs, for this small example min.edges is set to zero such that each tip clade corresponds to an individual locus
#' clusters <- extractClusters(ldna, min.edges=0, phi=0.25, rm.COCs=FALSE)
#' summary <- summaryLDna(ldna, clusters, LDmat)

summaryLDna <- function(ldna, clusters, LDmat){
  if(is.na(LDmat[2,1])) LDmat <- t(LDmat)
  p <- NULL
  for(i in 1:length(clusters)){
    fc <- names(clusters)[i]
    p <- c(p, as.vector(ldna$stats[, 2][ldna$stats[, 1] %in% fc]))
  }  
  Merge.at <- do.call('rbind', strsplit(p, "_"))[,2]
  
  nLoci <- as.vector(unlist(lapply(clusters, length)))
  
  #nE <- NULL
  #for(i in 1:length(clusters)){
  #  nE <- c(nE, unique(ldna$stats$nE[ldna$stats$cluster==names(clusters)[i]]))
  #}
  nE <- ldna$stats$nE[ldna$stats$cluster %in% names(clusters)]
  names(nE) <- ldna$stats$cluster[ldna$stats$cluster %in% names(clusters)]
  lambda <- ldna$stats$lambda[ldna$stats$cluster %in% names(clusters)]
  names(lambda) <- ldna$stats$cluster[ldna$stats$cluster %in% names(clusters)]

  std.err <- function(x) sd(x)/sqrt(length(x))
  Median.LD <- NULL
  MAD.LD <- NULL
  for(i in 1:length(names(clusters))){
    loci <- as.vector(unlist(clusters[i]))
    temp2 <- LDmat[which(rownames(LDmat) %in% loci), which(colnames(LDmat) %in% loci)]
    Median.LD <- c(Median.LD, signif(median(na.omit(as.vector(temp2))), digits=3))
    MAD.LD <- c(MAD.LD, signif(mad(as.vector(temp2), na.rm=TRUE, constant=1),3))
  }
  
  temp <- ldna$clusterfile[,colnames(ldna$clusterfile) %in% names(clusters)]
  nested <- matrix("SOC", ncol(temp), ncol(temp))
  for(i in 1:ncol(temp)){
    for(j in 1:ncol(temp)){
      if(i!=j & any(apply(cbind(temp[,i], temp[,j]),1, function(x) x[1]==TRUE & x[2]==TRUE))){
         nested[i,j] <- "COC"
      }
    }
  }
  nested[lower.tri(nested)] <- NA
  diag(nested) <- NA
  
  COCs <- as.vector(na.omit(colnames(temp)[apply(nested, 1, function(x) any(x=="COC"))]))
  SOCs <- colnames(temp)[!colnames(temp) %in% COCs]
  
  Type <- names(clusters)
  Type[which(Type %in% SOCs)] <- "SOC"
  Type[Type != "SOC"] <- "COC"
  names(Type) <- names(clusters)
  
  lambda[match(names(Type), names(lambda))]
  summary <- data.frame(Type, Merge.at, nLoci, nE=nE[match(names(Type), names(nE))], lambda=lambda[match(names(Type), names(lambda))], Median.LD, MAD.LD)
  return(summary)
}
