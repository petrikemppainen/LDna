#' Extracts LD clusters for linkage disequilibrium network analysis (LDna) 
#'
#' Identifies \emph{outlier clusters}, \emph{OCs}, and plots \code{"single linage clustering trees"} that describes cluster merging with decreasing LD threshold. Now also "branch traversal" is implemented.
#' 
#' If \code{plot.tree} and \code{plot.graph} are set to \code{TRUE}, \code{extractClusters} plots two graphs (default). The first shows all \eqn{\lambda}-values oredered from low to high and indicates which values are above \eqn{\lambda_{lim}}. Values corresponding to \emph{"selected outlier clusters", SOCs} are indicated in red and values corresponding to \emph{COCs} are indiced in blue. A \emph{COC} is defined as any \emph{OC} that contains loci from an \emph{OC} already extracted at a higher LD threshold. The second graph gives the tree illustrating cluster merger with decreasing LD threshold where nodes represent and node distance gives the LD threholds at which these events occur. Branches corresponding to \emph{SOCs} are indicated in red and those corresponding to \emph{COCs} are indiced in blue (if \code{rm.COCs=FALSE}).
#' 
#' Branch traversal means that as long as a branch contains outliers, the outlier clusters (SOC or COC) that is closest to the base of the branch is selected as the outlier cluster. This ensures that a wide range threshold values for \eqn{\lambda} gives similar results. 
#'
#' @param ldna Output from \code{\link{LDnaRaw}}
#' @param LDmat  Same LDmat as used for \code{\link{LDnaRaw}}
#' @param min.edges Minimum number of edges for a cluster that is shown as a branch in a tree.
#' @param phi Controls \eqn{\lambda_{lim}} which sets the threshold above which \eqn{\lambda} are considered as outliers. Default is two, values below this are not recommended.
#' @param rm.COCs If \code{TRUE} (default), automatically removes \emph{"compound oulier clusters" (COCs)}.
#' @param extract If \code{TRUE} (default), returns a list of cluster names. If \code{FALSE} only prints a tree with no \code{outlier clusters} indicated.
#' @param plot.tree If \code{TRUE} (default), plots tree. Has no effect if \code{extract=FALSE}.
#' @param plot.graph If \code{TRUE} (default), plots all \eqn{\lambda} ordered from low to high and indicates which values are above \eqn{\lambda_{lim}}.
#' @param branch.traversal If \code{TRUE} (default) "branch traversal" is implemented where clusters are extracted at the base of their branch or as far down the branch such that lambda is still above the required threshold.  
#' @param lambda.lim If not \code{NULL} gives a fixed value for \eqn{\lambda_{lim}}. Overrides any value passed by \code{phi}.
#' @keywords extractClusters
#' @seealso \code{\link{LDnaRaw}}, \code{\link{summaryLDna}} and \code{\link{plotLDnetwork}}
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, Christopher Knight \email{Chris.Knight@@manchester.ac.uk}
#' @return If extract=TRUE s list with two objects (1) a named list of vectors giving the loci for the extracted clusters and (2) a matrix indicating which 'single outlicer clusters' (\emph{SOCs}) are nested within which 'compound outlier clusters (\emph{COCs}). Else returns plots of the tree (if plot.tree=TRUE) and the distribution of \eqn{\lambda} values (if plot.graph=TRUE and extract=TRUE)
#' @examples
#' data(LDna)
#' 
#' ldna <- LDnaRaw(r2.baimaii_subs)
#' 
#' ## low lambda.lim without branch traversal, clusters extracted at high LD thresholds
#' par(mfcol=c(1,2))
#' clusters <- extractClusters(ldna, LDmat=r2.baimaii_subs, min.edges=15, lambda.lim = 0.1,extract=TRUE, plot.graph = TRUE, branch.traversal = FALSE)
#' ## low lambda.lim with branch traversal
#' # three different lambda.lim give exatly the same results after branch traversal; the "whole branch" is always selected. I recommend using branch traversal as it gives more consistent results across different data sets
#' par(mfcol=c(3,3))
#' clusters <- extractClusters(ldna, LDmat=r2.baimaii_subs, min.edges=15, lambda.lim = 0.1,extract=TRUE, plot.graph = FALSE, plot.tree=TRUE, branch.traversal = TRUE)
#' clusters <- extractClusters(ldna, LDmat=r2.baimaii_subs, min.edges=15, lambda.lim = 0.5,extract=TRUE, plot.graph = FALSE, plot.tree=TRUE, branch.traversal = TRUE)
#' clusters <- extractClusters(ldna, LDmat=r2.baimaii_subs, min.edges=15, lambda.lim = 1,extract=TRUE, plot.graph = FALSE, plot.tree=TRUE, branch.traversal = TRUE)
#' str(clusters)
#' # clusters are here 
#' clusters[[1]]
#' # matrix indicating nesting here 
#' clusters[[2]]
#' @export
extractClusters <- function(ldna, LDmat,min.edges=20, phi=2, lambda.lim=NULL, rm.COCs=TRUE, extract=TRUE, plot.tree=TRUE, plot.graph=TRUE, branch.traversal=TRUE){
  if(branch.traversal) rm.COCs=FALSE
  
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
    
    
    Ntips <- length(ldna$tree$tip.label)
    ## get outlier clusters
    clusters.out <- names(lambda_new)[which(lambda_new >= threshold)]
    if(identical(clusters.out, character(0))) stop("No outlier clusters, please decrease phi or lambda.lim")
    
    if(length(clusters.out)!=0){
      loci <- ldna$tree$tip.label
      # get SOCs and COCs
      temp <- ldna$clusterfile[,match(clusters.out, colnames(ldna$clusterfile))]
      if(is.matrix(temp)){
        nested <- matrix("SOC", ncol(temp), ncol(temp))
        for(i in 1:ncol(temp)){
          for(j in 1:ncol(temp)){
            if(i!=j & all(loci[temp[,j]] %in% loci[temp[,i]])){
              nested[i,j] <- "COC"
            }
          }
        }
        
        nested[lower.tri(nested)] <- NA
        COCs <- as.vector(na.omit(colnames(temp)[apply(nested, 1, function(x) any(x=="COC"))]))
        SOCs <- colnames(temp)[!colnames(temp) %in% COCs]
        diag(nested) <- NA
        
      }else{
        SOCs <- clusters.out
        COCs <- NA
        nested <- NULL
      }
      
    }else{
      
      print('No clusters to extract')
    } 
    
    dimnames(nested) <- list(clusters.out,clusters.out)
    
    
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
      plot(tree, show.tip.label=T, edge.width=3, edge.color=col, cex=1, tip.color=col.tip, root.edge=TRUE, underscore=T)
      axis(1, at=c(0,(1:10)*0.1))
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
    
    out <- list(clusters=loci, nested=nested)
    
    if(branch.traversal){
      
      out <- branch_traversal(ldna, clusters=out, Summary=as.data.table(summaryLDna(ldna=ldna,clusters = out,LDmat=LDmat)), min.edges = min.edges, lambda.lim = lambda.lim, plot.tree=plot.tree)
      return(out)
    }else{
      return(out)
    }
    
  }else{
    ifelse(max(LDmat, na.rm = TRUE)<1, 1,max(LDmat, na.rm = TRUE))
    plot(tree, edge.width=3,show.tip.label=F,edge.color="grey", cex=1,  root.edge=TRUE)
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
  tree <- ape::read.tree(text=newick)
}

branch_traversal <- function(ldna, clusters, Summary, min.edges = 20, lambda.lim = NULL,  plot.tree=FALSE, txt_tree1="Before traversal",txt_tree2="After traversal", cex=1){
  
  SOCs <- as.vector(Summary[Type=="SOC",Name])
  COCs <- as.vector(Summary[Type=="COC",Name])
  tree <- clusterPhylo(ldna, min.edges = min.edges)
  
  if(plot.tree){
    col <- rep("grey", length(tree$edge))
    tree$edge[tree$edge[,2] %in% which(tree$tip.label %in% SOCs),]
    distances <- tree$edge.length[tree$edge[,2] %in% which(tree$tip.label %in% SOCs)]
    clusters.temp <- tree$tip.label[tree$edge[,2][tree$edge[,2] %in% which(tree$tip.label %in% SOCs)]]
    keep.col <- clusters.temp[distances > 0]
    col[tree$edge[,2] %in% which(tree$tip.label %in% keep.col)] <- "red"
    col[tree$edge[,2] %in% tree$edge[,1][tree$edge[,2] %in% which(tree$tip.label %in% SOCs[!SOCs %in% keep.col])]] <- "red"
    col.tip <- rep("#00000000", length(tree$tip.label))
    
    col.tip[tree$tip.label %in% SOCs] <- "black"
    plot(tree, show.tip.label=T, edge.width=3, edge.color=col, cex=1, tip.color=col.tip, root.edge=TRUE, underscore=T)
    axis(1, at=c(0,(1:10)*0.1))
    if(!is.null(lambda.lim)){
      title(xlab="LD threshold", main=as.expression(bquote(lambda[lim]*plain("=")*.(lambda.lim)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
    }else{
      title(xlab="LD threshold", main=as.expression(bquote(varphi*plain("=")*.(phi)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
    }
    text(0,0,txt_tree1, cex=cex, adj = c(0,0))
    
  }
  
  n <- length(tree$tip.label)
  #cl <- SOCs[1]
  replace_cl <- do.call(rbind, lapply(SOCs, function(cl){
    # print(cl)
    anc <- cl 
    anc <- tree$edge[tree$edge[,2]==which(tree$tip.label==anc),1]
    keep <- prop.part(tree)[[anc - n]]
    clade <- drop.tip(tree, (1:n)[-keep], root.edge = 0, rooted = TRUE, collapse.singles = TRUE)
    # clade$tip.label
    go_on <- !any(SOCs %in% clade$tip.label[clade$tip.label != cl])
    branch <- NA
    while(go_on){
      # print(anc)
      anc <- tree$edge[tree$edge[,2]==anc,1]
      keep <- prop.part(tree)[[anc - n]]
      clade <- drop.tip(tree, (1:n)[-keep], root.edge = 0, rooted = TRUE, collapse.singles = TRUE)
      clade$tip.label
      go_on <- !any(SOCs %in% clade$tip.label[clade$tip.label != cl])
      if(go_on) branch <- unique(c(branch, COCs[COCs %in% clade$tip.label]))
      
    }
    # print(branch)
    
    if(!all(is.na(branch))){
      Min <- min(as.numeric(sapply(strsplit(branch, "_"), function(x) x[2])), na.rm = TRUE)
      
      if(length(na.omit(branch))!=0){
        new_soc <- branch[which(as.numeric(sapply(strsplit(branch, "_"), function(x) x[2]))==Min)]
        
        if(length(new_soc)>1){
          tmp <- lapply(new_soc, function(x){
            anc <- tree$edge[tree$edge[,2]==which(tree$tip.label==x),1]
            keep <- prop.part(tree)[[anc - n]]
            clade <- drop.tip(tree, (1:n)[-keep], root.edge = 0, rooted = TRUE, collapse.singles = TRUE)
            if(cl %in% clade$tip.label) return(c(cl,x[x!=cl]))  
          })
          if(length(which(sapply(tmp, length)!=0))>0) return(tmp[[which(sapply(tmp, length)!=0)]])
        }
        if(length(new_soc)==1){
          
          anc <- tree$edge[tree$edge[,2]==which(tree$tip.label==new_soc),1]
          keep <- prop.part(tree)[[anc - n]]
          clade <- drop.tip(tree, (1:n)[-keep], root.edge = 0, rooted = TRUE, collapse.singles = TRUE)
          if(cl %in% clade$tip.label) return(c(cl,new_soc[new_soc!=cl]))
          
        }
      }
    }
    
    
    
  }))
  
  
  SOCs[SOCs %in% replace_cl[,1]] <- replace_cl[,2]
  clusters.out <- SOCs
  
  if(length(clusters.out)!=0){
    loci <- ldna$tree$tip.label
    # get SOCs and COCs
    temp <- ldna$clusterfile[,match(clusters.out, colnames(ldna$clusterfile))]
    if(is.matrix(temp)){
      nested <- matrix("SOC", ncol(temp), ncol(temp))
      for(i in 1:ncol(temp)){
        for(j in 1:ncol(temp)){
          if(i!=j & all(loci[temp[,j]] %in% loci[temp[,i]])){
            nested[i,j] <- "COC"
          }
        }
      }
      
      nested[lower.tri(nested)] <- NA
      COCs <- as.vector(na.omit(colnames(temp)[apply(nested, 1, function(x) any(x=="COC"))]))
      SOCs <- colnames(temp)[!colnames(temp) %in% COCs]
      diag(nested) <- NA
      
    }else{
      SOCs <- clusters.out
      COCs <- NA
      nested <- NULL
    }
    
  }
  
  dimnames(nested) <- list(clusters.out,clusters.out)
  
  if(length(clusters.out)>1){
    out <- SOCs[order(-as.numeric(do.call('rbind', strsplit(SOCs, "_"))[,2]))]
    
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
  
  
  
  if(plot.tree){
    col <- rep("grey", length(tree$edge))
    tree$edge[tree$edge[,2] %in% which(tree$tip.label %in% SOCs),]
    distances <- tree$edge.length[tree$edge[,2] %in% which(tree$tip.label %in% SOCs)]
    clusters.temp <- tree$tip.label[tree$edge[,2][tree$edge[,2] %in% which(tree$tip.label %in% SOCs)]]
    keep.col <- clusters.temp[distances > 0]
    col[tree$edge[,2] %in% which(tree$tip.label %in% keep.col)] <- "red"
    col[tree$edge[,2] %in% tree$edge[,1][tree$edge[,2] %in% which(tree$tip.label %in% SOCs[!SOCs %in% keep.col])]] <- "red"
    col.tip <- rep("#00000000", length(tree$tip.label))
    
    col.tip[tree$tip.label %in% SOCs] <- "black"
    plot(tree, show.tip.label=T, edge.width=3, edge.color=col, cex=1, tip.color=col.tip, root.edge=TRUE, underscore=T)
    axis(1, at=c(0,(1:10)*0.1))
    if(!is.null(lambda.lim)){
      title(xlab="LD threshold", main=as.expression(bquote(lambda[lim]*plain("=")*.(lambda.lim)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
    }else{
      title(xlab="LD threshold", main=as.expression(bquote(varphi*plain("=")*.(phi)*"," ~~ "|E|"[min]*plain("=")* .(min.edges))))
    }
    text(0,0,txt_tree2, cex=cex, adj = c(0,0))
  }
  
  return(list(clusters=loci, nested=nested))
}