#' Extracts LD clusters for linkage disequilibrium network analysis (LDna) 
#'
#' Identifies \emph{outlier clusters}, \emph{OCs}, and plots \code{"single linage clustering trees"} that describes cluster merging with decreasing LD threshold.
#' 
#' If \code{plot.tree} and \code{plot.graph} are set to \code{TRUE}, \code{extractClusters} plots two graphs (default). The first shows all \eqn{\lambda}-values oredered from low to high and indicates which values are above \eqn{\lambda_{lim}}. Values corresponding to \emph{"selected outlier clusters", SOCs} are indicated in red and values corresponding to \emph{COCs} are indiced in blue. A \emph{COC} is defined as any \emph{OC} that contains loci from an \emph{OC} already extracted at a higher LD threshold. The second graph gives the tree illustrating cluster merger with decreasing LD threshold where nodes represent and node distance gives the LD threholds at which these events occur. Branches corresponding to \emph{SOCs} are indicated in red and those corresponding to \emph{COCs} are indiced in blue (if \code{rm.COCs=FALSE}).
#'
#' @param ldna Output from \code{\link{LDnaRaw}}
#' @param min.edges Minimum number of edges for a cluster that is shown as a branch in a tree.
#' @param phi Controls \eqn{\lambda_{lim}} which sets the threshold above which \eqn{\lambda} are considered as outliers. Default is two, values below this are not recommended.
#' @param rm.COCs If \code{TRUE} (default), automatically removes \emph{"compound oulier clusters" (COCs)}.
#' @param extract If \code{TRUE} (default), returns a list of cluster names. If \code{FALSE} only prints a tree with no \code{outlier clusters} indicated.
#' @param plot.tree If \code{TRUE} (default), plots tree. Has no effect if \code{extract=FALSE}.
#' @param plot.graph If \code{TRUE} (default), plots all \eqn{\lambda} ordered from low to high and indicates which values are above \eqn{\lambda_{lim}}.
#' @param lambda.lim If not \code{NULL} gives a fixed value for \eqn{\lambda_{lim}}. Overrides any value passed by \code{phi}.
#' @keywords extractClusters
#' @seealso \code{\link{LDnaRaw}}, \code{\link{summaryLDna}} and \code{\link{plotLDnetwork}}
#' @export
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, Christopher Knight \email{Chris.Knight@@manchester.ac.uk}
#' @return If extract=TRUE a named list of vectors giving the loci for the extracted clusters; otherwise returns plots of the tree (if plot.tree=TRUE) and the distribution of \eqn{\lambda} values (if plot.graph=TRUE and extract=TRUE)
#' @examples
#' # Simple upper diagonal LD matrix
#' LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.44, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.36, 0.24, NA, NA, NA, 0.48, 0.32, 0.2, 0.36, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.24, NA, NA, NA, NA, NA, NA, 0.44, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
#' # Calculate raw data
#' ldna <- LDnaRaw(LDmat)
#' # Extract clusters and plot graphs, for this example min.edges=0 such that each tip clade corresponds to an individual locus
#' par(mfcol=c(1,2))
#' clusters <- extractClusters(ldna, min.edges=0, phi=1)
#' clusters <- extractClusters(ldna, min.edges=0, phi=0.25, rm.COCs=FALSE)
#' # Larger data set
#' data(LDna)
#' ldna <- LDnaRaw(r2.baimaii_subs)
#' # Only print trees
#' clusters <- extractClusters(ldna, min.edges=15, extract=FALSE)
#' clusters <- extractClusters(ldna, min.edges=5, extract=FALSE)
#' # Different values of min.edges and phi can have crucial effects which clusters are extracted
#' clusters <- extractClusters(ldna, min.edges=15, phi=3)
#' clusters <- extractClusters(ldna, min.edges=10, phi=2)
#' # Include COCs
#' clusters <- extractClusters(ldna, min.edges=15, phi=5, rm.COCs=FALSE)
#' # Extract clusers without plottin graphs
#' clusters <- extractClusters(ldna, min.edges=15, phi=5, plot.tree=FALSE, plot.graph=FALSE)
#' # Set fixed value for lambda.lim
#' clusters <- extractClusters(ldna, min.edges=15, lambda.lim=2)

extractClusters <- function(ldna, min.edges=20, phi=2, lambda.lim=NULL, rm.COCs=TRUE, extract=TRUE, plot.tree=TRUE, plot.graph=TRUE){
  # Get file for tree and clusters above min.edges and their lambda values
  tree <- clusterPhylo(ldna, min.edges)
  if(extract){
    clusters <- tree$tip.label
    lambda_new <- ldna$stats$lambda[ldna$stats$nE >= min.edges][-1]
    names(lambda_new) <- as.vector(ldna$stats$cluster[ldna$stats$nE >= min.edges][-1])
    if(min.edges==0) lambda_new <- lambda_new[ldna$stats$nV!=1][-length(lambda_new[ldna$stats$nV!=1])]
    
    # Get thresholds. If lambda.lim is given, this takes precedence.
    if(!is.null(lambda.lim)){threshold <- lambda.lim
    }else{
      if(!is.null(phi)){
        threshold <- median(lambda_new)+mad(lambda_new, constant=phi)
      }  
    }
    
    #get outlier clusters
    clusters.out <- names(lambda_new)[which(lambda_new >= threshold)]
    if(identical(clusters.out, character(0))) stop("No outlier clusters, please decrease phi or lambda.lim")
    
    
    # get SOCs and COCs
    temp <- ldna$clusterfile[,colnames(ldna$clusterfile) %in% clusters.out]
    if(is.matrix(temp)){
      nested <- matrix(NA, ncol(temp), ncol(temp))
      for(i in 1:ncol(temp)){
        for(j in 1:ncol(temp)){
          if(i!=j & any(apply(cbind(temp[,i], temp[,j]),1, function(x) x[1]==TRUE & x[2]==TRUE))){
            nested[i,j] <- "COC"
          }
        }
      }
      nested[is.na(nested)] <- "SOC"
      nested[lower.tri(nested)] <- NA    
      COCs <- as.vector(na.omit(colnames(temp)[apply(nested, 1, function(x) any(x=="COC"))]))
      SOCs <- colnames(temp)[!colnames(temp) %in% COCs]
    }else{
      SOCs <- clusters.out
      COCs <- NA
    }
    
    if(plot.graph){
      col <- lambda_ord <- lambda_new[order(lambda_new)]
      #if(length(clusters.out)>1) 
      col[which(names(lambda_ord) %in% COCs)] <- "blue"
      col[which(names(lambda_ord) %in% SOCs)] <- "red"
      col[which(!col %in% c("red","blue"))] <- "black"
      plot(lambda_ord, col=col, ylab=expression(lambda))
      lines(c(1,length(lambda_ord)),c(threshold,threshold), col="red", lty=2) 
      text(length(lambda_ord)/3, y=threshold, as.expression(bquote(lambda[lim]*"="*.(signif(threshold,3)))), pos=3, adj=c(0,0))
      col.text <- lambda_ord
      if(rm.COCs==FALSE){
        if(length(clusters.out)>1) col.text[which(names(lambda_ord) %in% COCs)] <- "blue"
      } 
      col.text[which(names(lambda_ord) %in% SOCs)] <- "red"
      col.text[!col.text %in% c("blue","red")] <- "#00000000"
      text(lambda_ord, names(lambda_ord), pos=2, cex=0.75, col=col.text)
      if(!is.null(lambda.lim)){
        title(main=as.expression(bquote(lambda[lim]*plain("=")*.(lambda.lim)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
      }else{
        title(main=as.expression(bquote(varphi*plain("=")*.(phi)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
      }
    }
    
    # plot tree
    if(plot.tree){
      col <- rep("grey", length(tree$edge))
      if(rm.COCs==FALSE){
        distances <- tree$edge.length[tree$edge[,2] %in% which(tree$tip.label %in% clusters.out)]
        clusters.temp <- tree$tip.label[tree$edge[,2][tree$edge[,2] %in% which(tree$tip.label %in% clusters.out)]]
        keep.col <- clusters.temp[distances > 0]
        col[tree$edge[,2] %in% which(tree$tip.label %in% keep.col)] <- "blue"
        col[tree$edge[,2] %in% tree$edge[,1][tree$edge[,2] %in% which(tree$tip.label %in% clusters.out[!clusters.out %in% keep.col])]] <- "blue"
      }
      tree$edge[tree$edge[,2] %in% which(tree$tip.label %in% SOCs),]
      distances <- tree$edge.length[tree$edge[,2] %in% which(tree$tip.label %in% SOCs)]
      clusters.temp <- tree$tip.label[tree$edge[,2][tree$edge[,2] %in% which(tree$tip.label %in% SOCs)]]
      keep.col <- clusters.temp[distances > 0]
      col[tree$edge[,2] %in% which(tree$tip.label %in% keep.col)] <- "red"
      col[tree$edge[,2] %in% tree$edge[,1][tree$edge[,2] %in% which(tree$tip.label %in% SOCs[!SOCs %in% keep.col])]] <- "red"
      col.tip <- rep("#00000000", length(tree$tip.label))
      if(rm.COCs==FALSE){
        if(length(clusters.out)>1) col.tip[tree$tip.label %in% clusters.out] <- "blue"
      }
      col.tip[tree$tip.label %in% SOCs] <- "black"
      plot(tree, show.tip.label=T, edge.width=3, edge.color=col, cex=1, tip.color=col.tip, root.edge=TRUE, underscore=T,x.lim=1)
      #if(min.edges==0){
      axis(1, at=c(0,(1:10)*0.1))
      #}else{
      #  temp <- as.vector(ldna$stats[ldna$stats$nE > min.edges,2][-1])
      #  x <- round(10*max(as.numeric(do.call('rbind', strsplit(temp, "_", fixed=TRUE))[,2])),0)+1
      #  if(x>10) x <- 10
      #  axis(1, at=c(0,(1:x)*0.1))
      #}
      if(!is.null(lambda.lim)){
        title(xlab="LD threshold", main=as.expression(bquote(lambda[lim]*plain("=")*.(lambda.lim)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
      }else{
        title(xlab="LD threshold", main=as.expression(bquote(varphi*plain("=")*.(phi)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
      }
    }
    
    if(length(clusters.out)>1){
      
      if(rm.COCs==FALSE){out <- clusters.out[order(-as.numeric(do.call('rbind', strsplit(clusters.out, "_"))[,2]))]
      }else{out <- SOCs[order(-as.numeric(do.call('rbind', strsplit(SOCs, "_"))[,2]))]}
      
      temp <- ldna$clusterfile[,colnames(ldna$clusterfile) %in% out]
      if(is.matrix(temp)){
        temp <- temp[,order(as.numeric(do.call('rbind', strsplit(colnames(temp), "_", fixed=T))[,1]))]
        loci <- apply(temp, 2, function(x) rownames(temp)[x])          
      }else{
        loci <- list(names(temp)[temp])
      }
    }else{
      temp <- ldna$clusterfile[,colnames(ldna$clusterfile) %in% clusters.out]
      loci <- list(names(temp)[temp])
    }
    
    #rm(out)
    #rm(loci)
    return(loci)
    
  }else{
    plot(tree, edge.width=3,show.tip.label=F,edge.color="grey", cex=1,  root.edge=TRUE, x.lim=1)
    axis(1, at=c(0,(1:10)*0.1))
    title(xlab="LD threshold", main=as.expression(bquote("|E|"[min]*plain("=")* .(min.edges))))
  }
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
  tree <- read.tree(text=newick)
}