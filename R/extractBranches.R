#' Extracts LD clusters for linkage disequilibrium network analysis (LDna).
#'
#' Finds clusters of highly correlated loci (LD-clusters) by considering each branch in single-linkage clustering trees as a separate LD-cluster. The single-linkage can be plotted with the brances representing LD-clusters indicated.
#' 
#' In this implementation of defining LD-clusters, only the parameter \code{min.edges} needs to be specified. Low values of \code{min.edges} lead to "bushes" with many small "twigs" in the tree, each corresponding to a separate LD-cluster. These clusters contain fewer but more highly correlated sets of loci. Conversely, high values lead to fewer "thick" branches i.e. larger (in terms of the number of loci) but less strongly correlated LD-clusters.
#' \cr
#' \cr
#' The single-linkage clustering tree visualizing branches/LD-clusters can be printed by \code{plot.tree}=TRUE. 
#' \cr
#' \cr
#' With low values for \code{min.edges} an additional parameter \code{merge.min} can be used to group "branches" that merge above this threshold into a single LD-cluster; clusters that merge at high LD-threshold are so correlated, it makes no sense to regard them as separate LD-clusters.
#' 
#' This function works only with \code{\link{LDnaRaw2}} as shown in the examples
#' 
#' @param ldna Output from \code{\link{LDnaRaw}}
#' @param LDmat  Same LDmat as used for \code{\link{LDnaRaw}}
#' @param min.edges Minimum number of edges for a cluster that is shown as a branch in a tree.
#' @param merge.min Is the correlation threshold at which a clade is considered a separate LD-cluster, even if it contains more than two branches.
#' @param plot.tree If \code{TRUE} (default), plots tree.
#' @keywords extractBranches
#' @seealso \code{\link{LDnaRaw2}}
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, Christopher Knight \email{Chris.Knight@@manchester.ac.uk}
#' @return Returns a list containing a vectors of locus names belonging to a given cluster (branch)
#' @examples
#' data(LDna)
#' 
#' ldna <- LDnaRaw2(r2.baimaii_subs)
#' 
#' clusters <- extractBranches(ldna,min.edges=20,merge.min=0.8, plot.tree=TRUE) # default values
#' clusters
#' 
#' # with lower threshold for  \code{min.edges}
#' clusters <- extractBranches(ldna,min.edges=5,merge.min=0.8, plot.tree=TRUE) # default values
#' clusters
#' 
#' # with lower threshold for  \code{min.edges} and 
#' clusters <- extractBranches(ldna,min.edges=5,merge.min=0.2, plot.tree=TRUE) # default values
#' clusters
#' @export

extractBranches <- function(ldna, min.edges=20, merge.min=0.8, plot.tree=TRUE){
  
  
  # Get file for tree and clusters above min.edges and their lambda values
  tree <- clusterPhylo(ldna, min.edges)
  #plot(tree)
  
  clusterfile <- ldna$clusterfile
  
  
  stats <- data.table(ldna$stats)
  
  #suppressWarnings(clusterfile_red <- clusterfile[,!colnames(clusterfile) %chin% stats[merge_at_above >=merge.min & as.numeric(do.call(rbind, strsplit(stats$cluster, "_"))[,2])>merge.min, cluster]])
  clusterfile_red <- clusterfile[,colnames(clusterfile) %in% stats[nE>=min.edges,cluster] & colnames(clusterfile) %in% stats[merge_at_below<=merge.min,cluster]]
  
  SOCs <- clusters.out <- colnames(clusterfile_red)
  
  clusters.out <- unique(unlist(sapply(clusters.out, function(x) {
    anc <- stats[cluster %chin% x, parent_cluster]
    cand <- stats[parent_cluster %chin% anc, cluster]
    cand <- cand[!cand %chin% x]
    cand <- cand[cand %chin% clusters.out]
  })))
  
  if(length(clusters.out)!=0){
    loci <- ldna$tree$tip.label
    # get SOCs and COCs
    temp <- ldna$clusterfile[,colnames(ldna$clusterfile) %chin% clusters.out]
    if(is.matrix(temp)){
      nested <- matrix(NA, ncol(temp), ncol(temp))
      for(i in 1:ncol(temp)){
        for(j in 1:ncol(temp)){
          if(i!=j & all(loci[temp[,j]] %chin% loci[temp[,i]])){
            nested[i,j] <- "COC"
          }
        }
      }
      nested[is.na(nested)] <- "SOC"
      nested[lower.tri(nested)] <- NA
      COCs <- as.vector(na.omit(colnames(temp)[apply(nested, 1, function(x) any(x=="COC"))]))
      SOCs <- colnames(temp)[!colnames(temp) %chin% COCs]
      diag(nested) <- NA
      
    }else{
      SOCs <- clusters.out
      COCs <- NA
      nested <- NULL
    }
    
  }else{
    
    print('No clusters to extract')
  } 
  
  Ntips <- length(ldna$tree$tip.label)
  
  col <- rep("grey", length(tree$edge))
  
  
  
  if(plot.tree){
    distances <- tree$edge.length[tree$edge[,2] %in% which(tree$tip.label %in% SOCs)]
    clusters.temp <- tree$tip.label[tree$edge[,2][tree$edge[,2] %in% which(tree$tip.label %in% SOCs)]]
    keep.col <- clusters.temp[distances > 0]
    col[tree$edge[,2] %in% which(tree$tip.label %in% keep.col)] <- "red"
    col[tree$edge[,2] %in% tree$edge[,1][tree$edge[,2] %in% which(tree$tip.label %in% SOCs[!SOCs %in% keep.col])]] <- "red"
    col.tip <- rep("#00000000", length(tree$tip.label))
    
    col.tip[tree$tip.label %in% SOCs] <- "black"
    plot(tree, show.tip.label=TRUE, edge.width=3, edge.color=col, cex=0.5, tip.color=col.tip, root.edge=TRUE, underscore=TRUE)
    
    roof <- floor(ldna$stats[nE>min.edges,max(merge_at_above, na.rm = TRUE)]*10)
    axis(1, at=c(0,(1:roof)*0.1))
    text(0,1,as.expression(bquote("|E|"[min]*plain("=")* .(min.edges))),  adj = c(0,0))
  }
 
  clusters <- lapply(SOCs, function(x) names(which(clusterfile_red[,which(colnames(clusterfile_red)==x)])))
  names(clusters) <- SOCs
  return(clusters)
}
clusterPhylo <-  function(ldna, min.edges=0){
  d <- ldna$stats[ldna$stats$nE>=min.edges,]
  
  d$times<-rep(0,dim(d)[1])
  d$offspring<-as.list(rep(NA,dim(d)[1]))
  d$terminal<-rep(0,dim(d)[1])
  d$ancestor<-rep(0,dim(d)[1])
  
  row<-1
  rowOld<-0
  newick<-character()
  
  repeat{
    cluster<-as.character(d$cluster[row])
    d$offspring[[row]]<-which(as.character(d$parent_cluster)==cluster)
    offspNo<-sum(!is.na(d$offspring[[row]]))
    
    if(offspNo==0){d$terminal[row]<-1}
    
    if(d$terminal[row]==1){
      d$ancestor[row]<-rowOld
      if(d$ancestor[row]==0){
        newick<-paste("(",cluster,":0,:0);", sep="")
        break}
      newick<-paste(newick,cluster,":",d$distance[row], sep="")
      row<-d$ancestor[row]
      next}
    
    if(d$times[row]==0){
      newick<-paste(newick,"(", sep="")
      if(offspNo>1){newick<-paste(newick,"(", sep="")}
      d$ancestor[row]<-rowOld
      rowOld<-row
      d$times[row]<-d$times[row]+1
      row<-d$offspring[[row]][1]
      next}
    
    if(d$times[row]>0 & d$times[row]<offspNo){
      newick<-paste(newick,",", sep="")
      if((offspNo-d$times[row])>1){newick<-paste(newick,"(", sep="")}
      rowOld<-row
      d$times[row]<-d$times[row]+1
      row<-d$offspring[[row]][(d$times[row])]
      next}
    
    if(d$times[row]>0 & d$times[row]==offspNo){
      if(offspNo==1){
        newick<-paste(newick,",",cluster,":0):",d$distance[row], sep="")
      }else{
        newick<-paste(newick,paste(rep("):0",(offspNo-1)),sep="", collapse=""),",",cluster,":0):",d$distance[row], sep="", collapse="")            
      }
      if(d$ancestor[row]==0){
        newick<-paste(newick,";", sep="")
        break}
      row<-d$ancestor[row]
      next}
  }  
  tree <- ape::read.tree(text=newick)
}