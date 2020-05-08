#' Linkage disequilibrium (LD) network clustering, with pre-prepared edge lists as.
#'
#' Finds clusters of loci connected by high LD within non-overlapping windows and summerises this information using Principal Component Analysis, to be use e.g. in Genome Wide Association (GWA) analyses. For population genomic analyses it can also be used to filter non-infomative SNPs. Path to folder that contains edge lits (el's) calculated by other software needs to be specified. This implementetion is of course much faster since no LD values need to be calculated here.
#' 
#' Uses complete linkage clustering within non-overlapping windows of SNPs (within chromosomes) to find groups of correlated SNPs connected by high LD. 
#' This version el's of precalculated LD values as input. Assumes first column corresponds to locus 1 and the second column corresponds to locus 2 (of a pairwise comparison). The column that contain LD (r^2) values is specified with  \code{LD_column}
#' \cr
#' Window breakpoints are placed where r^2 between adjacent SNPs within a window (defined by \code{w1}) is below \code{LD_threshold1}. 
#' \code{nSNPs} determines the goal size for each window; smaller windows will produce faster compuation times. \cr
#' \cr
#' First, a graph object is created using all pairwise r^2 values between SNPs within each window, and edges (representing r^2 values) < LD_threshold1 are removed. This breaks the initial graph into sub-graphs. 
#' Removing edges like this and decomposing the graph is equivalent to a \code{sinlge linkage} clustering algorithm, with the draw back that clusters formed via single linkage clustering may be forced together due to single elements being close to each other, even though many of the elements in each cluster may be very distant to each other ('chaining-phenomenon').
#' For the (potentially long and elongated) sub-clusters produced in step one, a \code{complete linkage clustering} algorithim is performed (all edged need to be connected by a threshold value of r^2, resulting in more compact and highly connected clusters). \cr
#' \cr
#' With the \code{complete linkage clustering} trees clades with median LD>\code{LD_threshold2} are are recursively found. 
#' Thus \code{LD_threshold1} is typically set to lower than \code{LD_threshold1} and aims to break the initial network into sub-clusters, thus increasing compuational speed (i.e. it would also be possible to produce a complete linkage clustering tree for the whole window and then recursively find the relevant clades to extract). \cr
#' \cr
#' From each cluster as many Principal Components (\code{PCs}) as is required to explain >= \code{PC_threshold} are extracted. 
#' To increase compuational speed for large windows, or if a physical limit to how large clusters one is intrested in, window size (\code{w2}) for estimating r^2 can also be set to something less than the window size. This will in practice limit size of the clusters to (\code{w2}).
#' This algorithm is performed by chromosome and window, and can thus be parallelized (\code{cores}).\cr
#' \cr
#' It is also possible to plot LD-networks for different cluster solutions. For this you can specify a path for where the figure will be saved (\code{plot.network}). 
#' The parameter \code{threshold_net} determines the minimum value for each edge that is allowed in the network (i.e is equivalent to a single linkage clustering network). 
#' Each LD-cluster has a unique (arbitrary) color combination and black vertices reprsent loci that are not part of any cluster, i.e. should be considered independently in GWA analsyes. 
#' Minimum for \code{plot.network} is './', which exports figure to the working directory. Only works for data sets with up to 1000 SNPs and it is recommended that \code{threshold_net} is ~LD2.
#' 
#' @param EL_path Path to folder (only) containg relevant el's (one per chromosome/linkage group). For convenience the file name should be "name_of_chromosome"."el" e.g. chr1.el, 1.el, LG1.el, LG1_clean.el etc.
#' @param LD_column Index of the column that contain LD (r^2). The two first two columns must be locus 1 and 2 for each edge.
#' @param out_folder Path to folder where output is produced (default is './LDnCl_out/'). Proceed to use \code{Concat_files} to concatenate this data to a single file.
#' @param nSNPs Desired number of SNPs per window.
#' @param w1 Window size for defining putative recombination hotspots. Default is 10.
#' @param w2 Window size for estimating LD values (smaller is faster but may break up large clusters caused e.g. by inversion polymorphism).
#' @param LD_threshold1 Minimum LD value within a cluster
#' @param LD_threshold2 Minimum median LD within each cluster
#' @param PC_threshold Minimum cummulative amount of genetic variation explained by extracted PCs
#' @param cores Number of cores, default is 1
#' @param min.cl.size If you alrady know that you are only interested in loci that have LD>\code{LD_threshold2} with at least one other locus within the window specified by \code{w2}, \code{min.cl.size} can be set to 2. This greatly reduced the size of the output without greatly losing information e.g. for outlier analyeses, in particular in smaller data sets with otherwise high noise to signal ratios.
#' @param plot.network File name for plotting network. If \code{NULL} (default) no network is plotted.
#' @param threshold_net Threshold for edges when plotting network
#' @param to_do Vector with indexes for files in folder specified by \code{EL_path} that still need to be done. If \code{NULL} (default) all files will be processed
#' @keywords Linkage disequilibrium network clustering, complexity reduction
#' @seealso \code{\link{emmax_group}}
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, zitong.li \email{lizitong1985@@gmail.com}
#' @return Returns a list of five objects. These are saved as separate .rds files in the folder "LDnCl_out" (which is created if it does not already exist). If something goes wrong you can specify \code{to_do} with indexes for files in folder specified by \code{vcf_folder} that still need to be done. 
#' \code{cluster_summary} is a data frame that contains most of the relevant information for each cluster, inluding summary statistics (Median LD, MAD etc, see below). 
#' \code{cluster_PCs} includes a matrix individuals (rows) and PCs (columns) corresponding to each row in \code{cluster_summary}. This is what you use for subsequent GWA analyses. 
#' The third object (\code{clusters}) contains a list of locus names (chr:position; each entry corresponding to a row in \code{cluster_summary} for entries where PC=1). 
#' The fourth object is a vector of SNP positions for the highest intra cluster median LD ('maximally connected SNP', MCL) for each cluster, that can be used for downstream analyses
#' \cr 
#' \cr
#' Each .rds file contains only information for each chromosome but this data can be concatenated into a single file using \code{\link{Concat_files}}. 
#' \cr
#' The columns in file \code{cluster_summary} are:
#' \item{Chr}{Chromosome or linkage group identifer}
#' \item{Window}{Window identifier, recycled among chromosomes}
#' \item{Pos}{Mean position of SNPs in a cluster}
#' \item{Min}{Most downstream position of SNPs in a cluster}
#' \item{Max}{Most upstream position of SNPs in a cluster}
#' \item{Range}{Max-Min}
#' \item{nSNPs}{Number of SNPs in the cluster}
#' \item{Median}{Median pair wise LD between loci in a cluster}
#' \item{MAD}{Median absolute deviation among pair wise LD-values between loci in a cluster}
#' \item{PC}{PC identifier (first PC explaining most of the genetic variation)}
#' \item{PVE}{Cummulative proportion of variation explained by each PC}
#' \item{Grp}{Unique group identifier for PCs extracted from each cluster. Used by \code{\link{emmax_group}}}
#' @references Kemppainen, P., Knight, C. G., Sarma, D. K., Hlaing, T., Prakash, A., Maung Maung, Y. N., Walton, C. (2015). Linkage disequilibrium network analysis (LDna) gives a global view of chromosomal inversions, local adaptation and geographic structure. Molecular Ecology Resources, 15(5), 1031-1045. https://doi.org/10.1111/1755-0998.12369\cr
#' \cr
#' Li, Z., Kemppainen, P., Rastas, P., Merila, J. Linkage disequilibrium clustering-based approach for association mapping with tightly linked genome-wide data. Accepted to Molecular Ecology Resources.
#' @examples
#' ## example run, works basically the same as 'LDnClustering' but instead of defining 'snp' and 'map' files, you have to specify the path to a folder with el's. Works with only default values, if the correlct files are in "./EL_folder/"
#' \dontrun{
#' data(LDna)
#' ##create default folder
#' system("mkdir LD_EL") # create the default folder
#' ## put some edgelists in there
#' fwrite(ELs[[1]], "./LD_EL/chr1.el", sep = "\t")
#' fwrite(ELs[[2]], "./LD_EL/chr2.el", sep = "\t")
#' ## LD clustering using default parameters
#' LDnClusteringEL()
#' concat <- Concat_files("./LDnCl_out/") ## the default out folder
#' LDnCl_summary <- data.table(concat$cluster_summary)
#' LDnCl_summary
#' }
#' 
#' @export
#'

LDnClusteringEL <- function(EL_path = "./LD_EL", nSNPs=1000, w1=10, w2=100, LD_column = 4, out_folder="LDnCl_out", 
                               LD_threshold1 = 0.5, LD_threshold2 = 0.8, PC_threshold=0.8, plot.network=NULL, threshold_net=0.9, to_do = NULL, min.cl.size=2, cores=1){
  
  out <- list()
  if(is.null(to_do)){
    EL_files <- list.files(EL_path)  
  }else{
    EL_files <- list.files(EL_path)[to_do]
  }
  system(paste("mkdir ",out_folder), ignore.stderr = TRUE)
  
  
  for(i in 1:length(EL_files)){
    cat(paste0("Finding windows in file: ", EL_files[i], "\n"))
    
    
    chr <- strsplit(EL_files[i], ".", fixed = TRUE)[[1]][1]
    
    el <- fread(paste(EL_path, EL_files[i], sep = "/"))
    ## get ld between adjacent pairs of loci
    setnames(el, colnames(el)[c(1,2, LD_column)], c("from", "to", "r2"))
    
    el[,from := strip_white(from)]
    el[,to := strip_white(to)]

    el_adj <- as.matrix(el[!duplicated(from),.(from, to, r2)])
    snp.pos <- as.vector(unlist(el_adj[,1]))
    snp.pos <- c(snp.pos, as.vector(unlist(el_adj[nrow(el_adj),2])))
    
    nL <- length(snp.pos)
    nWindows <- as.integer(nL/nSNPs)+1
    
    temp <- which(slideFunct_max(el_adj[,3], w1, 1) < LD_threshold1)
    
    hotspots <- sapply(seq(1, nL, length.out = nWindows + 1), function(i) {which.min(abs(i - temp))})
    
    hotspots <- temp[hotspots]
    hotspots[length(hotspots)] <- nL
    hotspots[1] <- 0
    
    Windows <- lapply(1:(length(hotspots) - 1), function(x) {
      (1:nL)[(hotspots[x]+1):(hotspots[x + 1])]
    })
    cat(paste0("Number of windows: ", length(Windows), "; window  sizes: ", 
               paste(sapply(Windows, length), collapse = ":"), "\n"))
   
    out <- mclapply((1:length(Windows)), function(w){
      cat(paste0('Working on file ', EL_files[i], ', window ', w, '\n'))
      ## remove previous files, if present
      temp <- list.files()
      Window <- Windows[[w]]
      snp.pos_w <- snp.pos[Window]
      
      keep <- el$from %in% snp.pos_w & el$to %in% snp.pos_w
      g <- graph.edgelist(apply(el[keep,.(from, to)],2, function(o) as.character(o)), directed = FALSE)
      E(g)$weight <- as.numeric(el[keep,r2])
      
      g <- delete_edges(g, which(E(g)$weight<LD_threshold1))
    
      d_g <- decompose.graph(g)
      d_g <- d_g[sapply(d_g, function(x) length(V(x)$name)>1)]
      
     
      LDmat <- get.adjacency(g,attr = 'weight',names = TRUE,sparse = FALSE)
      diag(LDmat) <- NA
      LDmat <- t(LDmat)
      LDmat[upper.tri(LDmat)] <- NA
      
      PC_recursive <- function(clade){
        lapply(clade, function(x){
          
          if(x>ntips){
            cl <- ape::extract.clade(tree, x)$tip.label
          }else{
            cl <- tree$tip.label[x]
          }
          
          if(length(cl)>1){
            LDmat.part <- LDmat[which(snp.pos_w %in% cl), which(snp.pos_w %in% cl)]
            LDmat.part[upper.tri(LDmat.part)] <- t(LDmat.part)[upper.tri(LDmat.part)]
          
            Median.temp <- median(LDmat.part, na.rm  = TRUE)
            if(Median.temp > LD_threshold2 ){
              clusters.sub[[x]] <<- paste(chr, cl, sep=":")
              MCL.sub[[x]] <<- cl[which.max(apply(LDmat.part, 1, function(h) median(h, na.rm = TRUE)))]
              Median.sub[[x]] <<- Median.temp
              Mad.sub[[x]] <<- mad(LDmat.part, na.rm  = TRUE)
              
              
            }else{
              PC_recursive(tree$edge[,2][tree$edge[,1]==x])
            }
          }else{
            
            clusters.sub[[x]] <<- paste(chr, cl, sep=":")
            MCL.sub[[x]]  <<- cl
            Median.sub[[x]] <<- Mad.sub[[x]] <<- NA
          
          }
        })
      }
    
      clusters <- list()
      MCL <- list()
      Median <- list()
      Mad <- list()
      
      for(y in 1:length(d_g)){
       
        loci <- V(d_g[[y]])$name
        
        if(length(loci)>1){
          LDmat.part <- LDmat[which(snp.pos_w %in% loci), which(snp.pos_w %in% loci)]
          
          tree <- ape::as.phylo(hclust(as.dist(1-LDmat.part), method='complete'))
          tree$tip.label <- loci
          ntips <- length(tree$tip.label)
          
          clusters.sub <- list()
          MCL.sub <- list()
          Median.sub <- list()
          Mad.sub <- list()
          
          invisible(PC_recursive(tree$edge[,2][tree$edge[,1]==tree$edge[1,1]]))
          
          Keep <- !sapply(clusters.sub, is.null)
          clusters[[y]] <- clusters.sub[Keep]
          MCL[[y]] <- MCL.sub[Keep]
          Median[[y]] <- Median.sub[Keep]
          Mad[[y]] <- Mad.sub[Keep]
          
        }else{
          clusters[[y]] <- paste(chr, loci, sep=":")
          MCL[[y]]  <- loci
          Median[[y]] <- NA
          Mad[[y]] <- NA
        }
      }
      
     
      clusters <- flatten(clusters)
   
      MCL <- flatten(MCL)
      Mad <- unlist(flatten(Mad))
      Median <- unlist(flatten(Median))
      
      keep <- sapply(clusters, function(x) length(x))>=min.cl.size
      MCL <- MCL[keep]
      clusters <- clusters[keep]
      Mad <- Mad[keep]
      Median <- Median[keep]
      
      ## plot network if path for 'plot.network' and 'threshold_net' have been specified
      if(!is.null(plot.network)){
        ## plot network 
        cat('plotting network \n')
        
        nClust <- length(clusters)
        
        temp1 <- sample(rep(col_vector, ceiling(nClust/length((col_vector)))))
        
        
        Col <- unlist(lapply(1:length(clusters), function(x) {
          ifelse(sapply(clusters, length)[x]==1, return('black'),  return(rep(temp1[1:nClust][x], sapply(clusters, length)[x])))
        }))[match(snp.pos_w, unlist(clusters) )]
        
        
        temp2 <- sample(rep(col_vector[-1], ceiling(nClust/length((col_vector)))))
        
        frame.col <- unlist(lapply(1:length(clusters), function(x) {
          ifelse(sapply(clusters, length)[x]==1, return('black'),  return(rep(temp2[1:nClust][x], sapply(clusters, length)[x])))
        }))[match(colnames(LDmat), unlist(clusters) )]
        
        g <- graph.adjacency(LDmat, mode="lower", diag=FALSE, weighted=TRUE)
        
        V(g)$color <- Col
        V(g)$frame.color <- frame.col
        
        E(g)$weight <- round(E(g)$weight, 5)
        g <- delete.edges(g, which(E(g)$weight<=threshold_net))
        g <- delete.vertices(g, which(degree(g) == 0))
        
        par(mar=c(0,0,2,0))
        png(paste0(paste0(plot.network, 'Chr', i, '_window',w, '_LD2:', LD_threshold2), '.png'), res = 300, width = 6, height=6, units = 'in')
        plot.igraph(g, layout=layout.fruchterman.reingold, vertex.size=3, vertex.label.dist=NULL, edge.width=1, vertex.label=NA)
        title(main=paste(" @", threshold_net, sep="", '; Chr', i, '; window',w, '; LD2=', LD_threshold2))
        dev.off()
        
      }
      
      cat(paste0(length(MCL),' clusters and ', length(MCL), ' PCs extracted from of ', ncol(LDmat), ' original SNPs \n' ))
      
      return(list(clusters=clusters,  MCL=MCL, Median=Median,chr=chr))
    }, mc.cores=cores)
   
    cat('preparing output \n')
    clusters <- flatten(lapply(out, function(window){
      window$clusters
    }))
    
    MCL <- do.call(cbind, flatten(
      lapply(out, function(window){
        lapply(window$MCL, as.matrix)
      })
    ))
    
    Window <- unlist(lapply(1:length(out), function(window){
      rep(window, ncol(do.call(cbind, out[[window]][[2]])))
    }))
    
    Chr <- unlist(lapply(1:length(out), function(window){
      rep(out[[window]]$chr, ncol(do.call(cbind, out[[window]][[2]])))
    }))
    
    Median <- unlist(sapply(out, function(window){
      window$Median
    }))
    
    
    nSNPs2 <- sapply(clusters, function(x) length(x))
    Pos <- round(sapply(clusters, function(x) mean(as.numeric(sapply(strsplit(x, ":", fixed = TRUE), function(y) y[2])))))
    Min <- round(sapply(clusters, function(x) min(as.numeric(sapply(strsplit(x, ":", fixed = TRUE), function(y) y[2])))))
    Max <- round(sapply(clusters, function(x) max(as.numeric(sapply(strsplit(x, ":", fixed = TRUE), function(y) y[2])))))
    
    Range <- abs(Min-Max)
    
    cluster_summary <- data.frame(Chr, Window, Pos=Pos, Min, Max, Range, nSNPs=nSNPs2, Median=Median)
    
    colnames(cluster_summary) <- c('Chr', 'Window', 'Pos', 'Min', 'Max', 'Range', 'nSNPs', 'Median')
    saveRDS(list(cluster_summary=cluster_summary, clusters=clusters, MCL=MCL), file=paste0(out_folder,"/LDnCl:", chr, ".rds"))
    cat('Done \n')
  }
}  

Concat_files <- function(LDnCl_out){
  file_names <- list.files(LDnCl_out)
  
  LDnCl_combined <- list()
  LDnCl_combined$cluster_summary  <- do.call(rbind, lapply(file_names, function(x){
    readRDS(paste(LDnCl_out, x, sep = '/'))$cluster_summary
  }))
  
  LDnCl_combined$cluster_PCs <- do.call(cbind, lapply(file_names, function(x){
    readRDS(paste(LDnCl_out, x, sep = '/'))$cluster_PCs
  }))
  
  LDnCl_combined$clusters <- unlist(lapply(file_names, function(x){
    readRDS(paste(LDnCl_out, x, sep = '/'))$clusters
  }), recursive=FALSE)
  
  LDnCl_combined$MCL <- do.call(cbind, lapply(file_names, function(x){
    readRDS(paste(LDnCl_out, x, sep = '/'))$MCL
  }))
  
  LDnCl_combined$Cons <- do.call(cbind, lapply(file_names, function(x){
    readRDS(paste(LDnCl_out, x, sep = '/'))$Cons
  }))
  return(LDnCl_combined)
}
slideFunct_max <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=1)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- max(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}
