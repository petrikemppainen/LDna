#' Summerizes LD clusters
#'
#' Creates a data frame with summary information of LD clusters extracted by the function \code{\link{extractClusters}}
#'
#' @param ldna Output from \code{\link{LDnaRaw}}
#' @param clusters A file created by \code{\link{extractClusters}} that contains a list of all extracted clusters and their loci
#' @param LDmat Lower diagonal matrix of pairwise LD values
#' @keywords summaryLDna
#' @seealso \code{\link{extractClusters}}, \code{\link{LDnaRaw}}
#' @return Returns a data frame with each row corresponding to a cluster, in decreasing order with respect to highest LD threshold value at which they are present. Column \code{Type} specifies if a cluster is a "single outlier cluster" (SOC) or a "compound oulier cluster" (COC). Column "Merge.at" specfies the LD threshold at which the cluster merges. Column "Mean.LD" gives the average LD of all pairwise LD values between all loci in a cluster. The remaining columns are self exlanatory.
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @examples
#' # Simple upper diagonal LD matrix
#' LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.12, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.08, 0.4, NA, NA, NA, 0.48, 0.32, 0.04, 0.44, 0.36, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.32, NA, NA, NA, NA, NA, 0.72, 0.68, 0.28, NA, NA, NA, NA, NA, NA, 0.2, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
#' # Transform
#' # Get necessary data from  LD matrix
#' ldna <- LDnaRaw(LDmat)
#' # extract clusters and plot graphs, for this small example min.edges is set to zero such that each tip clade corresponds to an individual locus
#' clusters <- extractClusters(ldna, min.edges=0, phi=1, rm.COCs=FALSE)
#' summary <- summaryLDna(ldna, clusters, LDmat)

summaryLDna <- function(ldna, clusters, LDmat){
  
  p <- NULL
  for(i in 1:length(clusters)){
    fc <- names(clusters)[i]
    p <- c(p, as.vector(ldna$stats[, 2][ldna$stats[, 1] %in% fc]))
  }  
  Merge.at <- do.call('rbind', strsplit(p, "_"))[,2]
  
  nLoci <- as.vector(unlist(lapply(clusters, length)))
  
  #nE <- NULL
  #for(i in 1:length(clusters)){
  #  nE <- c(nE, unique(ldna$stats$nE[rownames(ldna$stats)==names(clusters)[i]]))
  #}
  nE <- ldna$stats$nE[rownames(ldna$stats) %in% names(clusters)]
  names(nE) <- rownames(ldna$stats)[rownames(ldna$stats) %in% names(clusters)]
  Lambda <- ldna$stats$Lambda[rownames(ldna$stats) %in% names(clusters)]
  names(Lambda) <- rownames(ldna$stats)[rownames(ldna$stats) %in% names(clusters)]

  std.err <- function(x) sd(x)/sqrt(length(x))
  Mean.LD <- NULL
  Mean.LD.SE <- NULL
  for(i in 1:length(names(clusters))){
    loci <- as.vector(unlist(clusters[i]))
    temp2 <- LDmat[which(rownames(LDmat) %in% loci), which(colnames(LDmat) %in% loci)]
    Mean.LD <- c(Mean.LD, signif(mean(na.omit(as.vector(temp2))), digits=3))
    Mean.LD.SE <- c(Mean.LD.SE, signif(std.err(na.omit(as.vector(temp2))),3))
  }
  
  temp <- ldna$clusterfile[,colnames(ldna$clusterfile) %in% names(clusters)]
  nested <- matrix("SOC", ncol(temp), ncol(temp))
  #nested <- rep("SOC", ncol(temp))
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
  
  Lambda[match(names(Type), names(Lambda))]
  summary <- data.frame(Type, Merge.at, nLoci, nE=nE[match(names(Type), names(nE))], Lambda=Lambda[match(names(Type), names(Lambda))], Mean.LD, Mean.LD.SE)
  return(summary)
}
