#' Plots linkage disequilibrium networks
#'
#' Allows for visual inspection of clusters identified by \code{\link{extractClusters}}.
#'
#' @keywords plotLDnetwork
#' @param ldna Output from \code{\link{LDnaRaw}}
#' @param el Edge list of all pairwise LD values for instance produced by \code{\link{mat2el}}
#' @param option \code{option=1} prints a full LD network from and edge list (\code{el}) and a specified LD threshold (\code{threshold}). \code{option=2} prints LD networks based on the files \code{clusters} and \code{summary} that contains information of extracted clusters.
#' @param threshold Specifies the LD threshold at which an LD network is plotted. Only required for option=1.
#' @param clusters Output from \code{\link{extractClusters}}.
#' @param summary Output from \code{\link{summaryLDna}}.
#' @param exl A list of locus names to be excluded from the LD networks (default is \code{NULL}).
#' @param full.network If \code{TRUE} (default), includes all loci in the LD networks, not recommended for large data sets.
#' @param include.parent If \code{full.network=FALSE} and \code{include.parent=TRUE} all loci from the parent cluster (after merger) are included. If \code{include.parent=FALSE} only the focal cluster is shown.
#' @param after.merger Whether to show LD networks at an LD threshold just before (\code{FALSE}) or just after (\code{TRUE}) merger.
#' @seealso \code{\link{LDnaRaw}}, \code{\link{extractClusters}}, \code{\link{mat2el}} and \code{\link{summaryLDna}}
#' @export
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @examples
#' ### Example with option=1
#' data(LDna)
#' plotLDnetwork(el=el, option=1, threshold=0.4)
#' ### Examples with option 2
#' par(mfcol=c(1,2))
#' ldna <- LDnaRaw(r2.baimaii_subs)
#' clusters <- extractClusters(ldna, min.edges=20, phi=5)
#' summary <- summaryLDna(ldna, clusters, r2.baimaii_subs)
#' el <- mat2el(r2.baimaii_subs)
#' #default settings with option=2
#' plotLDnetwork(ldna, el, option=2, clusters=clusters, summary=summary)
#' ## Other useful settings
#' # For large data sets
#' plotLDnetwork(ldna, el, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=FALSE, after.merger=FALSE)
#' # To visualise the merger
#' plotLDnetwork(ldna, el, option=2, clusters=clusters, summary=summary, full.network=TRUE, after.merger=FALSE)
#' # Or
#' plotLDnetwork(ldna, el, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=TRUE, after.merger=TRUE)
#' # To show that ususally several clusters are involved in most mergers
#' plotLDnetwork(ldna, el, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=TRUE, after.merger=FALSE)

plotLDnetwork <- function(ldna, el, option, threshold, clusters, summary,
exl=NULL, full.network=TRUE, include.parent=FALSE, after.merger=FALSE){
    if(option==1) g <- option1(el, threshold, exl);
    if(option==2) option2(ldna, el, clusters, summary, exl, full.network, include.parent, after.merger)
    if(option==1) return(g)
}

option1 <- function(el, threshold, exl){
    # exlude loci
    el <-el[!(el[,1] %in% exl),]
    el <-el[!(el[,2] %in% exl),]
    # only keep values above cut.off
    el.red <- el[which(as.numeric(el[,3])>threshold),]
    #calculate network
    g <- graph.edgelist(el.red[,1:2], directed=FALSE)
    plot.igraph(g, layout=layout.fruchterman.reingold, vertex.size=3, vertex.label.dist=NULL,
    vertex.color="grey", edge.width=1, vertex.label=NA, vertex.frame.color="grey")
    title(main=paste(" @", threshold, sep=""))
    return(g)
}

option2 <- function(ldna, el, clusters, summary, exl, full.network, include.parent, after.merger){
    col <- as.vector(summary$Type)
    col[col=="COC"] <- "blue"
    col[col=="SOC"] <- "red"
    for(i in 1:nrow(summary)){
        loci <- clusters[[i]]
        threshold <- as.numeric(as.vector(summary$Merge.at[rownames(summary) == names(clusters)[[i]]]))
        option2raw(ldna, el, exl, loci, col=col[i],
        full.network, threshold, include.parent, after.merger,
        names(clusters)[[i]])
    }
}

option2raw <- function(ldna, el, exl, loci, col, full.network, threshold, include.parent, after.merger, cluster.name){
    if(!is.null(exl)){
        el <- el[!(el[,1] %in% exl),]
        el <- el[!(el[,2] %in% exl),]
    }
    
    p <- as.vector(ldna$stats$parent_probe[ldna$stats$probe %in%  cluster.name])
    loci_p <- rownames(ldna$clusterfile)[ldna$clusterfile[,colnames(ldna$clusterfile) == p]]
    
    if(full.network==FALSE){
        if(include.parent==FALSE){
            el <-el[(el[,1] %in% loci),]
            el <-el[(el[,2] %in% loci),]
        }else{
            el <- el[(el[,1] %in% loci_p),]
            el <- el[(el[,2] %in% loci_p),]
        }
    }
    if(after.merger==TRUE) {
        p2 <- as.vector(ldna$stats$parent_probe[ldna$stats$probe %in%  p])
        threshold <- as.numeric(strsplit(p2, "_", fixed=T)[[1]][2])
    }
    el[,3]  <- round(as.numeric(el[,3]), 2)
    el.red <- el[which(as.numeric(el[,3])>threshold),]
    
    # only keep values above cut.off
    #calculate network
    g <- graph.edgelist(el.red[,1:2], directed=FALSE)
    
    #cluster color
    col.frame <- V(g)$name
    col.frame[which(col.frame %in% loci)] <- col
    col.frame[which(!col.frame %in% col)] <- "grey"
    
    plot.igraph(g, layout=layout.fruchterman.reingold, vertex.size=3, vertex.label.dist=NULL,
    vertex.color=col.frame, edge.width=1, vertex.label=NA, vertex.frame.color=col.frame)
    title(main=paste(cluster.name, " @", threshold, sep=""))
}