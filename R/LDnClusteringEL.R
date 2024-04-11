#' Linkage disequilibrium (LD) network clustering, starting from edge lists of pairwise LD-values as input.
#'
#' Finds clusters of loci connected by high LD within non-overlapping windows and returns a summary file, the resulting clusters and the names of the rSNPs (previously MCL), to be use e.g. in Genome Wide Association (GWA) or outlier analyses.
#' 
#' Uses single linkage clustering within non-overlapping windows of SNPs (within chromosomes) to find groups of correlated SNPs connected by high LD. This is done recursively within the single linkage clustering sub-trees (from root and up); as soon as a clade is reached where at least one locus is connected with all other loci within its clade above threshold \code{min_LD}, the algorithm stops. The default (0.7) produces clusters where the first PC typically explains >99 percent of the variation in each cluster.
#' 
#' Uses edge lists of LD values as input. The informative columns are specified by \code{LD_column}, where first and second specifies the columns for locus names and the third contains LD (r^2) values.
#' \cr
#' \cr
#' \code{nSNPs} determines the goal size for each window; smaller windows will produce faster computation times, but the smaller this size, the larger is the risk that you miss LD-clusters across breakpoints. However, if such a cluster is large, it will be split in two separate LD-clusters and they will be correlated in subsequent LDna steps.
#' \cr
#' \cr
#' If something goes wrong you can specify \code{to_do} with indexes for files in folder specified by \code{EL_path} that still need to be done. 
#' \cr
#' \cr
#' 
#' @param EL_path Path to folder (only) containing relevant el's (one per chromosome/linkage group). For convenience the file name should be "name_of_chromosome"."el" e.g. chr1.el, 1.el, LG1.el, LG1_clean.el etc. Note that it is unnecessary to use a window size of more than ca. 100 SNPs for the edge list (specified in whatever software you are using), otherwise you may have to wait for a long time...
#' @param nSNPs Desired number of SNPs per window.
#' @param columns Index of the column that contain LD (r^2). The two first two columns must be locus 1 and 2 for each edge.
#' @param out_folder Path to folder where output is produced (default is './LDnCl_out/'). Proceed to use \code{Concat_files} to concatenate this data to a single file.
#' @param min_LD Minimum LD value that at least one locus must be connected to all other loci in the recursive step.
#' @param cores Number of cores, default is 1
#' @param min.cl.size If 1 also singletons will be retained, which may be necessary for data sets with weak LD-structures (e.g. most loci are independent). Produced very large files though, so typically 2 is used to exclude all singleton clusters.
#' @param plot.network File name for plotting network. If \code{NULL} (default) no network is plotted.
#' @param threshold_net Threshold for edges when plotting network.
#' @param to_do Vector with indexes for files in folder specified by \code{EL_path} that still need to be done. If \code{NULL} (default) all files will be processed.
#' 
#' @keywords Linkage disequilibrium, LD, network analyses,complexity reduction
#' 
#' @seealso \code{\link{emmax_group}}
#' 
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, zitong.li \email{lizitong1985@@gmail.com}
#' 
#' @return Returns a list of three objects which are saved as ".rds" files in the folder specified by \code{out_folder}.
#' \code{cluster_summary} is a data frame that contains most of the relevant information for each cluster. 
#' \code{clusters} contains a list of locus names (chr:position; each entry corresponding to a row in \code{cluster_summary}.
#' \code{MCL} is a vector of names which which best represents the LD-cluster in downstream analyses ('maximally connected SNP', MCL aka rSNP).
#' \cr 
#' \cr
#' Each .rds file contains only information for each chromosome but they can be concatenated into a single file using \code{\link{Concat_files}}. 
#' \cr
#' \cr
#' The columns in file \code{cluster_summary} are:
#' Chr', 'Window', 'Pos', 'Min', 'Max', 'Range', 'nSNPs', 'Min_LD'
#' \item{Chr}{Chromosome or linkage group identifer}
#' \item{Window}{Window identifier, recycled among chromosomes}
#' \item{Pos}{Mean position of SNPs in a cluster}
#' \item{Min}{Most downstream position of SNPs in a cluster}
#' \item{Max}{Most upstream position of SNPs in a cluster}
#' \item{Range}{Max-Min}
#' \item{nSNPs}{Number of SNPs in the cluster}
#' \item{Min_LD}{the minimum LD between the rSNP/MCL and all other loci in its cluster}
#' @references Kemppainen, P., Knight, C. G., Sarma, D. K., Hlaing, T., Prakash, A., Maung Maung, Y. N., Walton, C. (2015). Linkage disequilibrium network analysis (LDna) gives a global view of chromosomal inversions, local adaptation and geographic structure. Molecular Ecology Resources, 15(5), 1031-1045. https://doi.org/10.1111/1755-0998.12369\cr
#' \cr
#' Li, Z., Kemppainen, P., Rastas, P., Merila, J. Linkage disequilibrium clustering-based approach for association mapping with tightly linked genome-wide data. Accepted to Molecular Ecology Resources.
#' @examples
#' \dontrun{
#' ## We will first create some example data to live in folder "LD_EL"
#' library(LDna)
#' data("LDna")
#' ## make directory for edge lists to live
#' system("mkdir LD_EL")
#' length(ELs) # edge lists for two chromosomes
#' ## write them in LD_EL folder
#' # the locus names need to be "Chr:Pos"
#' tmp <- as.data.table(ELs[[1]])
#' tmp[,V1:=paste("Chr1",V1,sep=":")]
#' tmp[,V2:=paste("Chr1",V2,sep=":")]
#' ELs[[1]] <- tmp
#' tmp <- as.data.table(ELs[[2]])
#' tmp[,V1:=paste("Chr2",V1,sep=":")]
#' tmp[,V2:=paste("Chr2",V2,sep=":")]
#' ELs[[2]] <- tmp
#' ## write the files to the EL folder
#' fwrite(ELs[[1]],file="LD_EL/Chr1.ld", row.names=FALSE,quote=FALSE)
#' fwrite(ELs[[2]],file="LD_EL/Chr2.ld", row.names=FALSE,quote=FALSE)
#' 
#' ## run LD-network clustering (LD-network complexity reduction)
#' LDnClusteringEL(EL_path = "./LD_EL/",cores = 10, min.cl.size = 2) ## no singleton clusters are kept
#' ## read in results
#' LDnC_res <- Concat_files("./LDnCl_out/")
#' cluster_summary<- as.data.table(LDnC_res$cluster_summary)
#' cluster_summary
#' cluster_summary[,hist(Min_LD)] ## distribution of minimum LD among any two loci within a cluster
#' cluster_summary[,table(nSNPs)]  ## distribution cluster sizes
#' cluster_summary[,plot(nSNPs,Min_LD)]  ## larger clusters tend to have lower minimum LD, those large clusters are from inversions
#' 
#' LDnC_res$clusters[cluster_summary[,which.max(nSNPs)]] ## this is the cluster with the most loci (e.g. putative inversion); the name is the MCL/rSNP
#' LDnC_res$MCL[cluster_summary[,which.max(nSNPs)]] ## the MCL/rSNP, i.e. the SNP that has the highest median LD with all other loci in this cluster 
#' ## and can be used to "represent" (hence "rSNP") this cluster in downstream analsyes.
#' ##  The alternative is to analyse the first PC as a forms of "synthetic multilocus genotypes"
#' }
#' @export

LDnClusteringEL <- function(EL_path = "./LD_EL/", nSNPs=1000,  columns = c(1,2,4), out_folder="LDnCl_out", 
                            min_LD = 0.7, plot.network=NULL, threshold_net=0.9, to_do = NULL, cores=1, min.cl.size=2){
  
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
    setnames(el, colnames(el)[columns], c("from", "to", "r2"))
    
    el[,from := strip_white(from)]
    el[,to := strip_white(to)]
    el[which(is.nan(r2)),r2:=0]
    el_adj <- as.matrix(el[!duplicated(from),.(from, to, r2)])
    snp.pos <- as.vector(unlist(el_adj[,1]))
    snp.pos <- c(snp.pos, as.vector(unlist(el_adj[nrow(el_adj),2])))
    
    nL <- length(snp.pos)
    nWindows <- as.integer(nL/nSNPs)+1
    
    temp <- which(slideFunct_max(as.numeric(el_adj[,3]), 1, 1) < 0.3)
    
    hotspots <- sapply(seq(1, nL, length.out = nWindows + 1), function(i) {which.min(abs(i - temp))})
    
    hotspots <- temp[hotspots]
    hotspots[length(hotspots)] <- nL
    hotspots[1] <- 0
    
    Windows <- lapply(1:(length(hotspots) - 1), function(x) {
      (1:nL)[(hotspots[x]+1):(hotspots[x + 1])]
    })
    cat(paste0("Number of windows: ", length(Windows), "; window  sizes: ", 
               paste(sapply(Windows, length), collapse = ":"), "\n"))
    
  
    #which.max(sapply(Windows,length))
    #w <- 48
    #cores <- 10
    out <- mclapply((1:length(Windows)), function(w){
      cat(paste0('Working on file ', EL_files[i], ', window ', w, '\n'))
      ## remove previous files, if present
      
      Window <- Windows[[w]]
      snp.pos_w <- snp.pos[Window]
      
      keep <- el$from %in% snp.pos_w & el$to %in% snp.pos_w
      g <- graph.edgelist(apply(el[keep,.(from, to)],2, function(o) as.character(o)), directed = FALSE)
      E(g)$weight <- as.numeric(el[keep,r2])
      
      g <- delete_edges(g, which(E(g)$weight<0.5))
      
      d_g <- decompose.graph(g)
      #d_g <- d_g[sapply(d_g, function(x) length(V(x)$name)>1)]
      
      
      LDmat <- el2mat_fast(el[keep,.(from, to, r2)])
      LDmat[upper.tri(LDmat)] <- t(LDmat)[upper.tri(LDmat)]
      
      PC_recursive <- function(clade){
        lapply(clade, function(x){
          
          if(x>ntips){
            cl <- ape::extract.clade(tree, x)$tip.label
          }else{
            cl <- tree$tip.label[x]
          }
          
          if(length(cl)>1){
            LDmat.part <- LDmat[which(locus_names %in% cl), which(locus_names %in% cl)]
            LDmat.part[LDmat.part==0] <- NA
            Min_LD <- apply(LDmat.part, 1, min, na.rm=TRUE)
            
            if(max(Min_LD) > min_LD){
              clusters.sub[[x]] <<- cl
              MCL.sub[[x]] <<- names(which.max(Min_LD))
              Median.sub[[x]] <<- max(Min_LD)
              
            }else{
              PC_recursive(tree$edge[,2][tree$edge[,1]==x])
            }
          }else{
            
            clusters.sub[[x]] <<- cl
            MCL.sub[[x]] <<-  cl
            Median.sub[[x]] <<- NA
            
          }
        })
      }
      
      clusters <- list()
      MCL <- list()
      Median <- list()
      
      locus_names <- colnames(LDmat)
      #y <- 1
      for(y in 1:length(d_g)){
        #print(y)
        loci <- V(d_g[[y]])$name
        
        LDmat.part <- LDmat[which(locus_names %in% loci), which(locus_names %in% loci)]
        if(length(loci)>1){
          
          # LDmat.part[1:10, 1:10]
          tree <- ape::as.phylo(hclust(as.dist(1-LDmat.part), method='single'))
          tree$tip.label <- loci
          ntips <- length(tree$tip.label)
          
          LDmat.part[LDmat.part==0] <- NA
          
          if(min(LDmat.part, na.rm = TRUE)<min_LD){
            clusters.sub <- list()
            MCL.sub <- list()
            Median.sub <- list()
            invisible(PC_recursive(tree$edge[,2][tree$edge[,1]==tree$edge[1,1]]))
            Keep <- !sapply(clusters.sub, is.null)
            clusters[[y]] <- as.vector(clusters.sub[Keep])
            MCL[[y]] <- MCL.sub[Keep]
            Median[[y]] <- Median.sub[Keep]
          }else{
            
            
            
            
            Median.temp <- median(LDmat.part, na.rm  = TRUE)
            #Min.temp <- min(LDmat.part, na.rm  = TRUE)
            
            clusters[[y]] <- as.vector(loci)
            MCL[[y]] <- loci[1]
            Median[[y]] <- Median.temp
            
          }
        }else{
          #Median.temp <- median(LDmat.part, na.rm  = TRUE)
          #LDmat.part <- LDmat[which(locus_names %in% loci), which(locus_names %in% loci)]
          clusters[[y]] <- loci
          MCL[[y]] <- loci
          Median[[y]] <- NA
          
        }
      }
      
      clusters <- flatten(clusters)
      MCL <- flatten(MCL)
      Median <- unlist(flatten(Median))
      
      keep <- which(sapply(clusters, function(x) length(x))>=min.cl.size)
      MCL <- MCL[keep]
      clusters <- clusters[keep]
      names(clusters) <- MCL
      Median <- Median[keep]
      
      cat(paste0(length(MCL),' clusters extracted from ', ncol(LDmat), ' original SNPs \n' ))
      
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
    
    colnames(cluster_summary) <- c('Chr', 'Window', 'Pos', 'Min', 'Max', 'Range', 'nSNPs', 'Min_LD')
    saveRDS(list(cluster_summary=cluster_summary, clusters=clusters, MCL=MCL), file=paste0(out_folder,"/LDnCl:", chr, ".rds"))
    cat('Done \n')
  }
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
strip_white <- function(x){trimws(x, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")}

el2mat_fast <- function(el){
  node_names <- strip_white(el[, unique(c(from,to))])
  n <- length(node_names)
  
  el2 <- cbind(chmatch(as.vector(el[,from]), node_names),
               chmatch(as.vector(el[,to]), node_names))
  
  mat<-matrix(0, n, n, dimnames = list(node_names, node_names))
  
  el2[which(apply(is.na(el2[,c(2,1)]), 1, any)),c(2,1)]
  
  mat[el2[,c(2,1)]] <- as.numeric(el[,r2])
  diag(mat) <- NA
  mat[upper.tri(mat)] <- NA
  mat
}


slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i]+window)])
  }
  return(result)
}

flatten <- function(x) {
  if (!inherits(x, "list")) return(list(x))
  else return(unlist(c(lapply(x, flatten)), recursive = FALSE))
}

col_vector <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", 
                "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", 
                "#66A61E", "#E6AB02", "#A6761D", "#666666", "#A6CEE3", "#1F78B4", 
                "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", 
                "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", 
                "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", 
                "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", 
                "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", 
                "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", 
                "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", 
                "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", 
                "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", 
                "#CCEBC5", "#FFED6F")
