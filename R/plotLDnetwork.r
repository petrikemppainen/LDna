#' Plots linkage disequilibrium networks
#'
#' Allows for visual inspection of clusters identified by \code{\link{extractClusters}}.
#'
#' See examples for more details
#'
#' @keywords plotLDnetwork
#' @param ldna Output from \code{\link{LDnaRaw}}
#' @param LDmat Lower diagonal LD matrix
#' @param option \code{option=1} prints a full LD network from and edge list (\code{LDmat}) and a specified LD threshold (\code{threshold}). \code{option=2} prints LD networks based on the files \code{clusters} and \code{summary} that contains information of extracted clusters.
#' @param threshold Specifies the LD threshold at which an LD network is plotted. Only required for option=1.
#' @param clusters Output from \code{\link{extractClusters}}.
#' @param summary Output from \code{\link{summaryLDna}}.
#' @param exl A list of locus names to be excluded from the LD networks (default is \code{NULL}).
#' @param full.network If \code{TRUE} (default), includes all loci in the LD networks, not recommended for large data sets.
#' @param include.parent If \code{full.network=FALSE} and \code{include.parent=TRUE} all loci from the parent cluster (after merger) are included. If \code{include.parent=FALSE} only the focal cluster is shown.
#' @param after.merger Whether to show LD networks at an LD threshold just before (\code{FALSE}) or just after (\code{TRUE}) merger.
#' @param graph.object Whether to output \code{graph object} when \code{option=1} (default=\code{FALSE}).
#' @seealso \code{\link{LDnaRaw}}, \code{\link{extractClusters}} and \code{\link{summaryLDna}}
#' @return If \code{option=1} and \code{graph.object=TRUE} the output is an igraph.object that can further be manipulated for custom networks (see \code{\link{igraph}} for details).
#' @export
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @examples
#' ### Example with option=1
#' data(LDna)
#' plotLDnetwork(LDmat=r2.baimaii_subs, option=1, threshold=0.4)
#' ### Examples with option 2
#' par(mfcol=c(1,2))
#' ldna <- LDnaRaw(r2.baimaii_subs)
#' clusters <- extractClusters(ldna, min.edges=20, phi=5)
#' summary <- summaryLDna(ldna, clusters, r2.baimaii_subs)
#' #default settings with option=2
#' plotLDnetwork(ldna, r2.baimaii_subs, option=2, clusters=clusters, summary=summary)
#' ## Other useful settings
#' # For large data sets
#' plotLDnetwork(ldna, r2.baimaii_subs, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=FALSE, after.merger=FALSE)
#' # To visualise the merger
#' plotLDnetwork(ldna, r2.baimaii_subs, option=2, clusters=clusters, summary=summary, full.network=TRUE, after.merger=TRUE)
#' # Or
#' plotLDnetwork(ldna, r2.baimaii_subs, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=TRUE, after.merger=TRUE)
#' # To show that ususally several clusters are involved in most mergers
#' plotLDnetwork(ldna, r2.baimaii_subs, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=TRUE, after.merger=FALSE)

plotLDnetwork <- function(ldna, LDmat, option, threshold, clusters, summary,
exl=NULL, full.network=TRUE, include.parent=FALSE, after.merger=FALSE, graph.object=FALSE){
    if(option==1) g <- option1(LDmat, threshold, exl);
    if(option==2) option2(ldna, LDmat, clusters, summary, exl, full.network, include.parent, after.merger)
    if(option==1 && graph.object) return(g)
}

option1 <- function(LDmat, threshold, exl){
    
    LDmat <- LDmat[!(rownames(LDmat) %in% exl),!(rownames(LDmat) %in% exl)]
    
    g <- graph.adjacency(LDmat, mode="lower", diag=FALSE, weighted=T)
    E(g)$weight <- round(E(g)$weight, 2)
    g <- delete.edges(g, which(E(g)$weight<=threshold))
    g <- delete.vertices(g, which(degree(g) == 0))
    
    plot.igraph(g, layout=layout.fruchterman.reingold, vertex.size=3, vertex.label.dist=NULL,
    vertex.color="grey", edge.width=1, vertex.label=NA, vertex.frame.color="grey")
    title(main=paste(" @", threshold, sep=""))
    return(g)
}

option2 <- function(ldna, LDmat, clusters, summary, exl, full.network, include.parent, after.merger){
    col <- as.vector(summary$Type)
    col[col=="COC"] <- "blue"
    col[col=="SOC"] <- "red"

    for(i in 1:nrow(summary)){
        threshold <- as.numeric(as.vector(summary$Merge.at[rownames(summary) == names(clusters)[[i]]]))
        option2raw(ldna, LDmat, exl, clusters[[i]], col=col[i],
        full.network, threshold, include.parent, after.merger,
        names(clusters)[[i]])
    }
}

option2raw <- function(ldna, LDmat, exl, loci, col, full.network, threshold, include.parent, after.merger, cluster.name){
    if(!is.null(exl)){
      LDmat <- LDmat[!(rownames(LDmat) %in% exl),!(rownames(LDmat) %in% exl)]
    }
    
    p <- as.vector(ldna$stats$parent_cluster[ldna$stats$cluster %in%  cluster.name])
    loci_p <- rownames(ldna$clusterfile)[ldna$clusterfile[,colnames(ldna$clusterfile) == p]]
    
    if(full.network==FALSE){
        if(include.parent==FALSE){
          LDmat <- LDmat[(rownames(LDmat) %in% loci), (rownames(LDmat) %in% loci)]
        }else{
          LDmat <- LDmat[(rownames(LDmat) %in% loci_p), (rownames(LDmat) %in% loci_p)]
        }
    }
    
    if(after.merger==TRUE) {
        p2 <- as.vector(ldna$stats$parent_cluster[ldna$stats$cluster %in%  p])
        if(p2=="root") {threshold <- threshold-0.01
        }else{threshold <- as.numeric(strsplit(p2, "_", fixed=T)[[1]][2])}
    }
    
    # only keep values above cut.off
    #calculate network
    g <- graph.adjacency(LDmat, mode="lower", diag=FALSE, weighted=T)
    E(g)$weight <- round(E(g)$weight, 2)
    g <- delete.edges(g, which(E(g)$weight<=threshold))
    g <- delete.vertices(g, which(degree(g) == 0))
    
    #cluster color
    col.frame <- V(g)$name
    col.frame[which(col.frame %in% loci)] <- col
    col.frame[which(!col.frame %in% col)] <- "grey"
    
    plot.igraph(g, layout=layout.fruchterman.reingold, vertex.size=3, vertex.label.dist=NULL,
    vertex.color=col.frame, edge.width=1, vertex.label=NA, vertex.frame.color=col.frame)
    title(main=paste(cluster.name, " @", threshold, sep=""))
}