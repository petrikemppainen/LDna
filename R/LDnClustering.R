#' Linkage disequilibrium (LD) network clustering
#'
#' Finds clusters of loci connected by high LD within non-overlapping windows and summerises this information using Principal Component Analysis, to be used in Genome Wide Association (GWA) analyses.
#' 
#' Uses complete linkage clustering within non-overlapping windows of SNPs to find clusters of SNPs connected by high LD. Window breakpoints are placed where r^2 between adjacent SNPs within a window (defined by \code{w1}) is below \code{LD_threshold1}. 
#'
#' @param snp A matrix with individuals as rows and bi-allelic SNP genotypes as columns. SNPs must be coded as 0,1 or 2 for the number of alleles of the reference nucleotide.
#' @param map A matrix with each row correspoinding a column in \code{snp}, column one corresponding to chromosome or linkage group and column two corresponding to physical position in the genome.
#' @param nSNPs Desired number of SNPs per window.
#' @param w1 Window size for defining putative recombination hotspots. Default is 10.
#' @param w2 Window size for estimating LD values (smaller is faster but may break up large clusters caused e.g. by inversion polymorphism).
#' @param LD_threshold1 Minimum LD value within a cluster
#' @param LD_threshold2 Minimum median LD within each cluster
#' @param PC_threshold Minimum cummulative amount of genetic variation explained by extracted PCs
#' @param verbose More detailed information of progress is printed on screen
#' @param plot.tree Plot complete linkage tree. 
#' @param mc.cores Number of cores for mclapply
#' @param plot.network File name for plotting network. If \code{NULL} (default) no network is plotted.
#' @param threshold_net Threshold for edges when plotting network.
#' @keywords Linkage disequilibrium network clustering, complexity reduction
#' @seealso \code{\link{emmax_group}}
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, zitong.li \email{lizitong1985@@gmail.com}
#' @return Returns a list of three objects. \code{map_cl} is a data frame that contains most of the relevant information for each cluster, inluding summary statistics (Median LD, MAD etc, see below). 
#' \code{PC_cl} includes a matrix of individuals (rows) and PCs (columns) corresponding to each row in map_pc. This is what you use for subsequent GWA analyses.
#' The third object (\code{cl_snps} contains a list of locus indexes (each entry corresponding to a row in map_pc for entries where PC=1).\cr
#' \cr
#' The columns in file \code{map_cl} are:
#' \item{Chr}{Chromosome or LG identifer}
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
#' \item{Grp}{Unique group identifier for PCs extracted from each cluster. Used by the emmax group algorithm.}
#' @references Kemppainen, P., Knight, C. G., Sarma, D. K., Hlaing, T., Prakash, A., Maung Maung, Y. N., … Walton, C. (2015). Linkage disequilibrium network analysis (LDna) gives a global view of chromosomal inversions, local adaptation and geographic structure. Molecular Ecology Resources, 15(5), 1031–1045. https://doi.org/10.1111/1755-0998.12369\cr
#' \cr
#' Li, Z., Kemppainen, P., Rastas, P., Merila, J. Linkage disequilibrium clustering-based approach for association mapping with tightly linked genome-wide data. Accepted to Molecular Ecology Resources.
#' @examples
#' ## Arabidopsis sample data (1000 first SNPs, Chr 1)
#' data(LDna)
#' 
#' t1 <- Sys.time()
#' data_0.1_0.3_old <- LDnClusteringOld(snp=snp_ara, 
#'                                      map=map_ara, 
#'                                      nSNPs = 1000, 
#'                                      w2 = 100, 
#'                                      LD_threshold1 = 0.1, 
#'                                      LD_threshold2 = 0.3, 
#'                                      verbose = TRUE)
#' t2 <- Sys.time()
#' data_0.1_0.3_new <- LDnClustering(snp=snp_ara, 
#'                                   map=map_ara, 
#'                                   nSNPs = 1000, 
#'                                   w2 = 100, 
#'                                   LD_threshold1 = 0.1, 
#'                                   LD_threshold2 = 0.3, 
#'                                   verbose = TRUE)
#' t3 <- Sys.time()
#' 
#' difftime(t3, t2)[[1]]/difftime(t2, t1)[[1]]
#' 
#' t1 <- Sys.time()
#' data_0.3_0.5_old <- LDnClusteringOld(snp=snp_ara, 
#'                                      map=map_ara, 
#'                                      nSNPs = 1000, 
#'                                      w2 = 100, 
#'                                      LD_threshold1 = 0.3, 
#'                                      LD_threshold2 = 0.5, 
#'                                      verbose = TRUE)
#' t2 <- Sys.time()
#' data_0.3_0.5_new <- LDnClustering(snp=snp_ara, 
#'                                   map=map_ara, 
#'                                   nSNPs = 1000, 
#'                                   w2 = 100, 
#'                                   LD_threshold1 = 0.3, 
#'                                   LD_threshold2 = 0.5, 
#'                                   verbose = TRUE)
#' t3 <- Sys.time()
#' difftime(t3, t2)[[1]]/difftime(t2, t1)[[1]]
#' 

LDnClustering <- function(snp, map, nSNPs=1000, w1=10, w2=100 ,LD_threshold1 = 0.1, LD_threshold2 = 0.3, PC_threshold=0.8, verbose=TRUE, plot.tree=FALSE, mc.cores=1, plot.network=NULL, threshold_net=0.9){
  
  out <- list()
  
  snp.id <- 1:nrow(map)
  
  snp.pos.org <- map[,2]
  #for each chromosome
  for(i in unique(map[,1])){
    
    cat("Finding windows \n")
    Chr <- map[, 1] == i
    snp.pos <- snp.pos.org[Chr]
    snp.id.chr <- snp.id[Chr]
    nL <- length(snp.pos)
    
    ## use snpGDS to estiamte r^2
    gds_name <- paste0(gsub(" ", "", Sys.time()), ".gds")
    snpgdsCreateGeno(gds_name, genmat = snp[, Chr], sample.id = 1:nrow(snp), snp.pos = snp.id.chr, snpfirstdim = FALSE)
    file <- snpgdsOpen(gds_name)
    LDmat <- snpgdsLDMat(file, method = "r", slide = 1, verbose = FALSE)$LD^2
    snpgdsClose(file)
    system(paste0("rm ", gds_name))
    
    
    ## find hotspots and windows
    nWindows <- as.integer(nL/nSNPs)
    
    
    if (nWindows > 1) {
      temp <- which(slideFunct(LDmat[1, ], w1) < LD_threshold1)
      
      hotspots <- sapply(seq(1, nL, length.out = nWindows + 1), function(i) {which.min(abs(i - temp))})
      
      hotspots <- temp[hotspots]
      hotspots[length(hotspots)] <- nL
      hotspots[1] <- 0
      
      Windows <- lapply(1:(length(hotspots) - 1), function(x) {
        snp.id.chr[(hotspots[x]+1):(hotspots[x + 1])]
      })
      
    }else {
      Windows <- list(snp.id.chr)
    }
    
    cat(paste0("Number of windows: ", length(Windows), "; window  sizes: ", 
               paste(sapply(Windows, length), collapse = ":"), "\n"))
    
    ## for each window
    #w <- 1
    out[[i]] <- mclapply(1:length(Windows), function(w){
      cat(paste('Working on chromosome', i ,', window', w, '\n'))
      
      Window <- Windows[[w]]
      snp.id.bin <- snp.id[which(snp.id%in%Window)]
      snp_bin  <- snp[,which(snp.id%in%Window)]
      
      
      snpgdsCreateGeno(paste0(w,gds_name), genmat = snp_bin, sample.id = 1: nrow(snp) , snp.id = snp.id.bin, snpfirstdim=FALSE) ### sample ID should be 1: how many individuals you have
      file <- snpgdsOpen(paste0(w,gds_name))
      
      # estimate LD within a sliding window
      if(w2>length(snp.id.bin)) w2 = -1 ## if window size longer that number of snps, estiamte LD for all pairwise comparisons
      
      MAT <- snpgdsLDMat(file, method="r",slide=w2, verbose=FALSE)
      snpgdsClose(file)
      system(paste0('rm ', paste0(w,gds_name)))
      
      
      if(w2!=-1){
        LDmat <- data.table(remove=rep(NA, length(MAT$snp.id)))
        for(x in 1:(length(MAT$snp.id))){
          LDmat[(1:w2+x)[!is.nan(MAT$LD[,x])],as.character(x) := MAT$LD[!is.nan(MAT$LD[,x]),x]^2]
        }
        LDmat <- as.matrix(LDmat)
        LDmat <- LDmat[,-1]
      }else{
        for(x in 1:(length(MAT$snp.id))){
          LDmat <- MAT$LD^2
        }
      }
      
      ## preparare LDmat
      LDmat[is.na(LDmat)]<- 0
      LDmat[upper.tri(LDmat)] <- NA
      diag(LDmat) <- NA
      colnames(LDmat) <- snp.id.bin
      rownames(LDmat) <- snp.id.bin
      
      ## pre-clustering of network; remove edges below LD_threshold1
      g <- graph.adjacency(LDmat, mode="lower", diag=FALSE, weighted=TRUE)
      g <- delete_edges(g, which(E(g)$weight<LD_threshold1))
      d_g <- decompose.graph(g)
      
      ## the recursive function; finds the clades where median LD is above LD_threshold2
      
      PC_recursive <- function(clade){
        lapply(clade, function(x){
          
          if(x>ntips){
            cl <- extract.clade(tree, x)$tip.label
          }else{
            cl <- tree$tip.label[x]
          }
          
          if(length(cl)>1){
            LDmat.part <- LDmat[which(rownames(LDmat) %in% cl), which(colnames(LDmat) %in% cl)]
            LDmat.part[upper.tri(LDmat.part)] <- t(LDmat.part[lower.tri(LDmat.part)])
            
            if(median(LDmat.part, na.rm  = TRUE) > LD_threshold2){
              clusters.sub[[x]] <<- cl
              PC_sc <- PC_score(snp_bin[,snp.id.bin %in% cl], PC_threshold)
              PVE.sub[[x]] <<- PC_sc[[2]]
              PCs.sub[[x]] <<- PC_sc[[1]]
            }else{
              PC_recursive(tree$edge[,2][tree$edge[,1]==x])
            }
          }else{
            clusters.sub[[x]] <<- cl
            PVE.sub[[x]] <<- 1
            PCs.sub[[x]] <<- as.matrix(snp_bin[,snp.id.bin %in% cl])
          }
        })
      }
      
      clusters <- list()
      PVE <- list()
      PCs <- list()
      for(y in 1:length(d_g)){
        loci <- V(d_g[[y]])$name
        if(length(loci)>1){
          LDmat.part <- LDmat[loci, loci]
          tree <- as.phylo(hclust(as.dist(1-LDmat.part), method='complete'))
          tree$tip.label <- loci
          
          ntips <- length(tree$tip.label)
          clusters.sub <- list()
          PVE.sub <- list()
          PCs.sub <- list()
          
          invisible(PC_recursive(tree$edge[1,1]))
          clusters[[y]] <- clusters.sub[!sapply(clusters.sub, is.null)]
          PVE[[y]] <- PVE.sub[!sapply(PVE.sub, is.null)]
          PCs[[y]] <- PCs.sub[!sapply(PCs.sub, is.null)]
          
        }else{
          clusters[[y]] <- loci
          PVE[[y]] <- 1
          PCs[[y]] <- as.matrix(snp_bin[,snp.id.bin %in% loci])
          
        }
      }
      
      clusters <- flatten(clusters)
      PVE <- flatten(PVE)
      PCs <- flatten(PCs)
      
      ## get Median and MAD
      stats <- do.call(rbind, lapply(clusters, function(x){
        loci <- as.vector(unlist(x))
        temp <- LDmat[which(rownames(LDmat) %in% loci), which(colnames(LDmat) %in% loci)]
        c(Median=signif(median(na.omit(as.vector(temp))), digits=5), MAD=signif(mad(as.vector(temp), na.rm=TRUE, constant=1),3))
      }))
      
      ## plot network if path for 'plot.network' and 'threshold_net' have been specified
      if(!is.null(plot.network)){
        ## plot network 
        cat('plotting network \n')
        
        clusters.temp <- clusters[sapply(clusters, length) != 1]
        nClust <- length(clusters)
        
        temp1 <- sample(rep(col_vector, ceiling(nClust/length((col_vector)))))
        
        Col <- unlist(lapply(1:length(clusters.temp), function(x) {
          ifelse(sapply(clusters.temp, length)[x]==1, return('black'),  return(rep(temp1[1:nClust][x], sapply(clusters.temp, length)[x])))
        }))[match(colnames(LDmat), unlist(clusters.temp) )]
        
        
        temp2 <- sample(rep(col_vector[-1], ceiling(nClust/length((col_vector)))))
        
        frame.col <- unlist(lapply(1:length(clusters.temp), function(x) {
          ifelse(sapply(clusters.temp, length)[x]==1, return('black'),  return(rep(temp2[1:nClust][x], sapply(clusters.temp, length)[x])))
        }))[match(colnames(LDmat), unlist(clusters.temp) )]
        
        par(mar=c(0,0,2,0))
        png(paste0(paste(plot.network, i, b, sep='_'), '.png'), res = 300, width = 6, height=6, units = 'in')
        plotLDnetwork(ldna,LDmat,option=1, threshold=threshold_net, col=Col, frame.color = frame.col)
        dev.off()
      }
      
      cat(paste0(length(PCs),' clusters and ', sum(sapply(PCs, ncol)), ' PCs extracted from of ', ncol(LDmat), ' original SNPs, on average explaining ', signif(mean(sapply(PVE, max)), 2), '% of the variation \n' ))
      
      return(list(clusters=clusters, PCs=PCs, PVE=PVE, stats=stats))
      
    }, mc.cores=mc.cores)
  }
  
  
  ### prepare output
  
  cat('preparing output \n')
  if(length(out)==1)out <- c(list(NULL), out)
  out <- out[sapply(out, length)>0]
  
  clusters <- flatten(lapply(out, function(chromosome){
    lapply(chromosome, function(window){
      window[[1]]
    })
  }))
  
  clusters <- lapply(clusters, as.numeric)
  
  cluster_PCs <- flatten(lapply(out, function(chromosome){
    lapply(chromosome, function(window){
      window[[2]]
    })
  }))
  
  Chr <- unlist(lapply(1:length(out), function(chromosome){
    lapply(1:length(out[[chromosome]]), function(window){
      rep(unique(map[,1])[chromosome], ncol(do.call(cbind, out[[chromosome]][[window]][[2]])))
    })
  }))
  
  Window <- unlist(lapply(1:length(out), function(chromosome){
    lapply(1:length(out[[chromosome]]), function(window){
      rep(window, ncol(do.call(cbind, out[[chromosome]][[window]][[2]])))
    })
  }))
  
  Pos <- unlist(sapply(out, function(chromosome){
    sapply(chromosome, function(window){
      rep(sapply(window[[1]], function(x) mean(map[as.numeric(x),2])), sapply(window[[2]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)))
    })
  }))
  
  PC <- unlist(sapply(out, function(chromosome){
    sapply(chromosome, function(window){
      sapply(sapply(window[[2]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)), function(bin) 1:bin)
    })
  }))
  
  nSNPs <- unlist(sapply(out, function(chromosome){
    sapply(chromosome, function(window){
      rep(sapply(window[[1]], length), sapply(window[[2]], function(bin) ncol(bin)))
    })
  }))
  
  PVE <- unlist(sapply(out, function(chromosome){
    sapply(chromosome, function(window){
      rep(window[[3]], sapply(window[[2]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)))
    })
  }))
  
  PVE <- unlist(lapply(out, function(chromosome){
    lapply(chromosome, function(window){
      window[[3]]
    })
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
  
  Median <-  unlist(sapply(out, function(chromosome){
    sapply(chromosome, function(window){
      rep(window[[4]][,1], sapply(window[[2]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)))
    })
  }))
  
  MAD <-  unlist(sapply(out, function(chromosome){
    sapply(chromosome, function(window){
      rep(window[[4]][,2], sapply(window[[2]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)))
    })
  }))
  
  Min <-  sapply(clusters, function(x) min(map[x,2]))[Grp]
  Max <-  sapply(clusters, function(x) max(map[x,2]))[Grp]
  Range <- abs(Min-Max)
  cluster_summary <- data.frame(data.table(Chr, Window, Pos=Pos, Min, Max, Range, nSNPs=nSNPs, Median=Median, MAD=MAD, PC=PC, PVE, Grp))
  colnames(cluster_summary) <- c('Chr', 'Window', 'Pos', 'Min', 'Max', 'Range', 'nSNPs', 'Median', 'MAD', 'PC', 'PVE', 'Grp')
  PC_data <- do.call(cbind, cluster_PCs)
  
  return(list(cluster_summary=cluster_summary, cluster_PCs=cluster_PCs, clusters=clusters))
}

PC_score <- function(x_A, PC_threshold){
  cMean_A <- colMeans(x_A)
  for (i in 1:dim(x_A)[2]) x_A[, i] <- x_A[, i] - cMean_A[i]
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
