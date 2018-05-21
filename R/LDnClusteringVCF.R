#' Linkage disequilibrium (LD) network clustering, with vcf files as input
#'
#' Finds clusters of loci connected by high LD within non-overlapping windows and summerises this information using Principal Component Analysis, to be use in Genome Wide Association (GWA) analyses. Takes .vcf files as input.
#' 
#' Uses complete linkage clustering within non-overlapping windows of SNPs (within chromosomes) to find groups of correlated SNPs connected by high LD. 
#' This version takes .vcf files as input and uses \code{vcftools} to estimate LD. Assumes that vcftools can be called from command line as 'vcftools', thus make sure that it has the correct permission (use chmod 755 if necessary) \cr
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
#' This algorithm is performed by chromosome and window, and can thus be parallelized (\code{mc.cores}).\cr
#' \cr
#' It is also possible to plot LD-networks for different cluster solutions. For this you can specify a path for where the figure will be saved (\code{plot.network}). 
#' The parameter \code{threshold_net} determines the minimum value for each edge that is allowed in the network (i.e is equivalent to a single linkage clustering network). 
#' Each LD-cluster has a unique (arbitrary) color combination and black vertices reprsent loci that are not part of any cluster, i.e. should be considered independently in GWA analsyes. 
#' Minimum for \code{plot.network} is './', which exports figure to the working directory. Only works for data sets with up to 1000 SNPs and it is recommended that \code{threshold_net} is ~LD2.
#' 
#' @param vcf_folder Path to folder (only) containg relevant vcf files (one per chromosome).
#' @param ld_measure Either '---geno-r2' or '---hap-r2', depending if you have phased data or not.
#' @param maf Minor allele frequency; SNPs need to be bi-allellic.
#' @param nSNPs Desired number of SNPs per window.
#' @param w1 Window size for defining putative recombination hotspots. Default is 10.
#' @param w2 Window size for estimating LD values (smaller is faster but may break up large clusters caused e.g. by inversion polymorphism).
#' @param LD_threshold1 Minimum LD value within a cluster
#' @param LD_threshold2 Minimum median LD within each cluster
#' @param PC_threshold Minimum cummulative amount of genetic variation explained by extracted PCs
#' @param mc.cores Number of cores for mclapply
#' @param plot.network File name for plotting network. If \code{NULL} (default) no network is plotted.
#' @param threshold_net Threshold for edges when plotting network
#' @param to_do Vector with indexes for files in folder specified by \code{vcf_folder} that still need to be done
#' @keywords Linkage disequilibrium network clustering, complexity reduction
#' @seealso \code{\link{emmax_group}}
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, zitong.li \email{lizitong1985@@gmail.com}
#' @return Returns a list of five objects. These are saved as separate .rds files in the folder "LDnCl_out" (which is created if it does not already exist). If something goes wrong you can specify \code{to_do} with indexes for files in folder specified by \code{vcf_folder} that still need to be done. 
#' \code{cluster_summary} is a data frame that contains most of the relevant information for each cluster, inluding summary statistics (Median LD, MAD etc, see below). 
#' \code{cluster_PCs} includes a matrix individuals (rows) and PCs (columns) corresponding to each row in \code{cluster_summary}. This is what you use for subsequent GWA analyses. 
#' The third object (\code{clusters}) contains a list of locus names (chr:position; each entry corresponding to a row in \code{cluster_summary} for entries where PC=1). 
#' The fourth and fifth objects are files matrixes similar to \code{cluster_PCs} except there will be one column per cluster with either the SNP with the highest intra cluster median LD ('maximally connected SNP', MCL) or the consensus SNP (Cons), respectively. 
#' See example code for details. These last two files are expected to work better when LD thresholds are higher and SNPs in clusters or more strongly correlated, but no simualtions have been done to test this yet. 
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
#' ## example run, works basically the same as 'LDnClustering' but instead of defining 'snp' and 'map' files, you have to specify the path to a folder with '.vcf' files.
#' \dontrun{
#' data_0.3_0.5 <- LDnClustering(vcf_folder="./vcf", 
#'                                   ld_measure='--geno-r2',
#'                                   nSNPs = 1000, 
#'                                   w2 = 100, 
#'                                   LD_threshold1 = 0.3, 
#'                                   LD_threshold2 = 0.5,
#'                                   PC_threshold = 0.8)
#' }
#'                   
#' @export
#' 
LDnClusteringVCF <- LDnClusteringVCF <- function(vcf_folder = "./vcf", ld_measure = '--geno-r2', nSNPs=1000, w1=10, w2=100 ,
                                                 LD_threshold1 = 0.5, LD_threshold2 = 0.7, PC_threshold=0.8, maf=0.05, mc.cores=1, plot.network=NULL, threshold_net=0.9, to_do = to_do){
  
  out <- list()
  vcf_files <- list.files(vcf_folder)[to_do]
  temp <- list.files()
  if(!any(grepl('LDnCl_out', temp))) system("mkdir LDnCl_out")
  
  for(i in 1:length(vcf_files)){
    cat(paste0("Finding windows in file: ", vcf_files[i], "\n"))
    
    temp <- list.files()
    if(any(grepl('ldchr', temp))) system(paste0("rm ", paste0("./", temp[grep('ldchr', temp)])))
    
    
    ## get bi-allelic loci
    system(paste0("vcftools --vcf ", vcf_folder, vcf_files[[i]], " ", "--maf ", maf, " --max-alleles 2 --min-alleles 2 --recode --recode-INFO-all"), ignore.stderr = TRUE)
    system(paste0("vcftools --vcf  out.recode.vcf --ld-window 1 ", ld_measure, " --out ldchr"), ignore.stderr = TRUE)
    
    ## get ld between adjacent pairs of loci
    
    ## get LD within chromosome to determine windows
    temp <- list.files()
    ld <- fread(paste0("./", temp[grep('ldchr', temp)]))
    system(paste0("rm ", paste0("./", temp[grep('ldchr', temp)])))
    
    chr <- unlist(ld[1,1])
    snp.pos <- as.vector(unlist(ld[,2]))
    snp.pos <- c(snp.pos, as.vector(unlist(ld[nrow(ld),3])))
    
    nL <- length(snp.pos)
    LDmat <- ld$`R^2`
    nWindows <- as.integer(nL/nSNPs)+1
    
    temp <- which(slideFunct(LDmat, w1) < LD_threshold1)
    hotspots <- sapply(seq(1, nL, length.out = nWindows + 1), function(i) {which.min(abs(i - temp))})
    
    hotspots <- temp[hotspots]
    hotspots[length(hotspots)] <- nL
    hotspots[1] <- 0
    
    Windows <- lapply(1:(length(hotspots) - 1), function(x) {
      (1:nL)[(hotspots[x]+1):(hotspots[x + 1])]
    })
    cat(paste0("Number of windows: ", length(Windows), "; window  sizes: ", 
               paste(sapply(Windows, length), collapse = ":"), "\n"))
    
    out <- lapply((1:length(Windows)), function(w){
      cat(paste0('Working on file ', vcf_files[i], ', window ', w, '\n'))
      ## remove previous files, if present
      temp <- list.files()
      if(any(grepl('ldw', temp))) system(paste0("rm ", paste0("./", temp[grep('ldw', temp)])))
      
      Window <- Windows[[w]]
      
      system(paste0("vcftools --vcf  out.recode.vcf ", " --from-bp ",  snp.pos[Window[1]], " --to-bp ", snp.pos[Window[length(Window)]], " --chr ", chr," --out subs --recode --recode-INFO-all"), ignore.stderr = TRUE)
      system(paste0("vcftools --vcf  subs.recode.vcf --ld-window ", w2, "  --out ldw ",ld_measure), ignore.stderr = TRUE)
      system("vcftools --vcf subs.recode.vcf --012 --out data", ignore.stderr = TRUE)
      
      system("sed -i -e 's/-1/NA/g' data.012")
      snp_w <- as.matrix(fread('data.012',header=FALSE))[,-1]
      snp.pos_w <-  snp.pos[Window]
      
      temp <- list.files()
      if(any(grepl('ldw.log', temp))) system(paste0("rm ", paste0("./", temp[grep('ldw.log', temp)])))
      ldw <- fread(paste0("./", temp[grep('ldw', temp)]))
      
      ldw$`R^2`[which(is.na(ldw$`R^2`))] <- 0
      
      system(paste0("vcftools --vcf  out.recode.vcf ", " --from-bp ",  snp.pos[Window[1]], " --to-bp ", snp.pos[Window[length(Window)]], " --chr ", chr," --out subs --recode --recode-INFO-all"), ignore.stderr = TRUE)
      
      system(paste0("rm ", paste0("./", temp[grep('ldw', temp)])))
      
      el <- as.matrix(ldw[,.(POS1, POS2, `R^2`)])
      g <- graph.edgelist(apply(el[,1:2],2, function(o) as.character(as.integer(o))), directed = FALSE)
      E(g)$weight <- as.numeric(el[,3])
      
      g <- delete_edges(g, which(E(g)$weight<LD_threshold1))
      d_g <- decompose.graph(g)
      
      
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
            if(Median.temp > LD_threshold2){
              snp_cl <- snp_w[,which(snp.pos_w %in% cl)]
              
              clusters.sub[[x]] <<- paste(chr, cl, sep=":")
              PC_sc <- PC_score(snp_w[,snp.pos_w %in% cl], PC_threshold)
              PVE.sub[[x]] <<- PC_sc[[2]]
              PCs.sub[[x]] <<- PC_sc[[1]]
              MCL.sub[[x]] <<- snp_cl[,which.max(apply(LDmat.part, 1, function(h) median(h, na.rm = TRUE)))]
              Median.sub[[x]] <<- Median.temp
              Mad.sub[[x]] <<- mad(LDmat.part, na.rm  = TRUE)
              # polarise sequences
              temp <- do.call(cbind,lapply(2:ncol(snp_cl), function(h){
                if(!all(snp_cl[,1]==snp_cl[,h])){
                  if(cor.test(snp_cl[,1], snp_cl[,h])$estimate<0){
                    return(ifelse(snp_cl[,h]==2, 0,2))
                  }else{
                    return(rev(snp_cl[,h]))
                  }
                }else{
                  return(rev(snp_cl[,h]))
                }
                
              }))
              
              # get consensus
              Cons.sub[[x]] <<- phyclust::find.consensus(t(temp))
            }else{
              PC_recursive(tree$edge[,2][tree$edge[,1]==x])
            }
          }else{
            clusters.sub[[x]] <<- paste(chr, cl, sep=":")
            PVE.sub[[x]] <<- 1
            Cons.sub[[x]] <<- MCL.sub[[x]] <<- PCs.sub[[x]] <<- as.matrix(snp_w[,snp.pos_w %in% cl])
            Median.sub[[x]] <<- Mad.sub[[x]] <<- NA
          }
        })
      }
      
      
      clusters <- list()
      PVE <- list()
      PCs <- list()
      MCL <- list()
      Cons <- list()
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
          PVE.sub <- list()
          PCs.sub <- list()
          MCL.sub <- list()
          Cons.sub <- list()
          Median.sub <- list()
          Mad.sub <- list()
          
          invisible(PC_recursive(tree$edge[,2][tree$edge[,1]==tree$edge[1,1]]))
          Keep <- !sapply(clusters.sub, is.null)
          clusters[[y]] <- clusters.sub[Keep]
          PVE[[y]] <- PVE.sub[Keep]
          PCs[[y]] <- PCs.sub[Keep]
          MCL[[y]] <- MCL.sub[Keep]
          Cons[[y]] <- Cons.sub[Keep]
          Median[[y]] <- Median.sub[Keep]
          Mad[[y]] <- Mad.sub[Keep]
          
        }else{
          clusters[[y]] <- paste(chr, loci, sep=":")
          PVE[[y]] <- 1
          Cons[[y]] <- MCL[[y]] <- PCs[[y]] <- as.matrix(snp_w[,snp.pos_w %in% loci])
          Median[[y]] <- NA
          Mad[[y]] <- NA
          
          
        }
      }
      
      
      clusters <- flatten(clusters)
      PVE <- flatten(PVE)
      PCs <- flatten(PCs)
      MCL <- flatten(MCL)
      Cons <- flatten(Cons)
      Mad <- unlist(flatten(Mad))
      Median <- unlist(flatten(Median))
      
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
      
      cat(paste0(length(PCs),' clusters and ', sum(sapply(PCs, ncol)), ' PCs extracted from of ', ncol(LDmat), ' original SNPs, on average explaining ', signif(mean(sapply(PVE, max)), 2)*100, '% of the variation \n' ))
      
      return(list(clusters=clusters, PCs=PCs, PVE=PVE, MCL=MCL, Cons=Cons, Median=Median, Mad=Mad, chr=chr))
    })
    
    cat('preparing output \n')
    
    length(out)
    clusters <- flatten(lapply(out, function(window){
      window[[1]]
    }))
    
    MCL <- do.call(cbind, flatten(
      lapply(out, function(window){
        lapply(window[[4]], as.matrix)
      })
    ))
    
    
    Cons <-  do.call(cbind, flatten(
      lapply(out, function(window){
        lapply(window[[5]], as.matrix)
      })
    ))
    
    
    cluster_PCs <- do.call(cbind, flatten(lapply(out, function(window){
      window[[2]]
    })))
    
    
    Window <- unlist(lapply(1:length(out), function(window){
      rep(window, ncol(do.call(cbind, out[[window]][[2]])))
    }))
    
    
    PC <- unlist(sapply(out, function(window){
      sapply(sapply(window[[2]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)), function(bin) 1:bin)
    }))
    
    PVE <- unlist(lapply(out, function(window){
      window[[3]]
    }))
    
    Grp <- rep(NA, length(PC))
    j <- 0
    for(x in 1:length(PC)){
      if(PC[x]==1){
        j <- j + 1
        Grp[x] <- j
      }else{
        Grp[x] <- j
      }
    }
    
    Chr <- unlist(lapply(1:length(out), function(window){
      rep(out[[window]][[8]], ncol(do.call(cbind, out[[window]][[2]])))
    }))
    
    Median <- unlist(sapply(out, function(window){
      window[[6]]
    }))[Grp]
    
    MAD <- unlist(sapply(out, function(window){
      window[[7]]
    }))[Grp]
    
    nSNPs2 <- sapply(clusters, function(x) length(x))[Grp]
    Pos <- round(sapply(clusters, function(x) mean(as.numeric(sapply(strsplit(x, ":", fixed = TRUE), function(y) y[2])))))[Grp]
    Min <- round(sapply(clusters, function(x) min(as.numeric(sapply(strsplit(x, ":", fixed = TRUE), function(y) y[2])))))[Grp]
    Max <- round(sapply(clusters, function(x) max(as.numeric(sapply(strsplit(x, ":", fixed = TRUE), function(y) y[2])))))[Grp]
    
    Range <- abs(Min-Max)
    
    cluster_summary <- data.frame(Chr, Window, Pos=Pos, Min, Max, Range, nSNPs=nSNPs2, Median=Median, MAD=MAD, PC=PC, PVE, Grp)
    
    colnames(cluster_summary) <- c('Chr', 'Window', 'Pos', 'Min', 'Max', 'Range', 'nSNPs', 'Median', 'MAD', 'PC', 'PVE', 'Grp')
    saveRDS(list(cluster_summary=cluster_summary, cluster_PCs=cluster_PCs, clusters=clusters, MCL=MCL, Cons=Cons), file=paste0("./LDnCl_out/LDnCl:", chr, ".rds"))
    cat('Done \n')
  }
  
}  

PC_score <- function(x_A, PC_threshold){
  cMean_A <- colMeans(x_A, na.rm = TRUE)
  for (i in 1:dim(x_A)[2]) x_A[, i] <- x_A[, i] - cMean_A[i]
  if(any(is.na(x_A))){
    invisible(lapply(1:ncol(x_A), function(k){
      x_A[is.na(x_A[,k]),k] <<- mean(x_A[,k], na.rm = TRUE)
    }))
  }
  eigen_A <- eigen(var(x_A))
  scores_A <- x_A %*% eigen_A$vectors
  percentage_A = (cumsum(eigen_A$values))/sum(eigen_A$values)
  a = percentage_A < PC_threshold
  a[which(!a)[1]] <- TRUE
  zeta = as.matrix(scores_A[, a])
  return(list(zeta, percentage_A[a]))
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


