#' Linkage disequilibrium (LD) network binning
#'
#' Finds clusters of loci connected by high LD within a sliding window summerizes genetic variation using principal component analyses
#' 
#' If \code{plot.tree} and \code{plot.graph} are set to \code{TRUE}, \code{extractClusters} plots two graphs (default). The first shows all \eqn{\lambda}-values oredered from low to high and indicates which values are above \eqn{\lambda_{lim}}. Values corresponding to \emph{"selected outlier clusters", SOCs} are indicated in red and values corresponding to \emph{COCs} are indiced in blue. A \emph{COC} is defined as any \emph{OC} that contains loci from an \emph{OC} already extracted at a higher LD threshold. The second graph gives the tree illustrating cluster merger with decreasing LD threshold where nodes represent and node distance gives the LD threholds at which these events occur. Branches corresponding to \emph{SOCs} are indicated in red and those corresponding to \emph{COCs} are indiced in blue (if \code{rm.COCs=FALSE}).
#'
#' @param snp A matrix with individuals as rows and bi-allelic SNP genotypes as columns. SNPs must be coded as 0,1 or 2 for the number of alleles of the reference nucleotide.
#' @param mapsnp A matrix with each row correspoinding a column in \code{snp}, column one corresponding to chromosome or linkage group and column two corresponding to physical position in the genome.
#' @param nSNPs Desired number of SNPs per bin.
#' @param w1 Window size for defining putative recombination hotspots.
#' @param w2 Window size for estimating LD values.
#' @param LD_threshold1 Minimum LD value within a cluster
#' @param LD_threshold2 Minimum median LD within each cluster, and for the second step of clustering
#' @param PC_threshold Minimum cummulative amount of genetic variation explained by extracted PCs
#' @param verbose More detailed information of progress is outputted
#' @param plot.tree Plot tree. Default \code{FALSE}
#' @param mc.cores number of cores for mclapply
#' @param plot.network File name for plotting network. If \code{NULL} (default) no network is plotted.
#' @param threshold_net Threshold for edges when plotting network.
#' @keywords LDnB
#' @seealso \code{\link{LDnaRaw}}, \code{\link{summaryLDna}} and \code{\link{extractClusters}}
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, zitong.li
#' @return Returns a list of three objects. \code{mapbin} contains map information for each PC extracted from the data, PCbin=PCbin, bin_snps=bin_snps \cr
#' test
#' @examples
#' # To come

LDnBin <- function(snp, mapsnp, nSNPs=1000, w1=10, w2=100 ,LD_threshold1 = 0.1, LD_threshold2 = 0.3, PC_threshold=0.8, verbose=FALSE, plot.tree=FALSE, mc.cores=1, plot.network=NULL, threshold_net=0.9){
  
  out <- list()
  #for each chromosome
  # i <- 12
  
  snp.id <- 1:nrow(mapsnp)
  
  snp.pos.org <- mapsnp[,2]
  
  for(i in unique(mapsnp[,1])){
    cat('Finding windows \n')
    snp.pos <- snp.pos.org[mapsnp[,1]==i]
    
    temp1 <- 10^max(unlist(lapply(strsplit(as.character(snp.pos), '.', fixed = TRUE), nchar)))
    
    snp.pos <- unlist(lapply(unique(snp.pos), function(x){
      temp <- snp.pos[snp.pos==x]
      
      
      temp+signif(((1:length(temp)/length(temp))-1/length(temp))/temp1
                  , log10(temp1)+2)
    }))
    
    gds_name <- paste0(gsub(' ', '', Sys.time()), '.gds')
    
    snpgdsCreateGeno(gds_name, genmat = snp[,mapsnp[,1]==i] , sample.id = 1: nrow(snp) , snp.pos = snp.id[mapsnp[,1]==i], snpfirstdim=FALSE) ### sample ID should be 1: how many individuals you have
    file <- snpgdsOpen(gds_name)
    
    ## get hotspots
    LDmat <- snpgdsLDMat(file, method="r",slide=1, verbose=FALSE)$LD^2
    
    snpgdsClose(file)
    system(paste0('rm ', gds_name))
    
    nWindows <- as.integer(length(snp.pos)/nSNPs)
    
    #plot(slideFunct(LDmat[1,], w1))
    
    
    if(nWindows>1){
      temp <- which(slideFunct(LDmat[1,], w1)<LD_threshold1)
      
      hotspots <- sapply(seq(1, length(snp.pos), length.out = nWindows+1)[-1], function(i){
        which.min(abs(i-temp))
      })
      hotspots <- c(0, temp[hotspots])
      
      Windows <- lapply(1:(length(hotspots)-1), function(x){
        snp.id[mapsnp[,1]==i][(hotspots[x]+1):(hotspots[x+1])]
      })
    }else{
      Windows <- list(snp.id[mapsnp[,1]==i])
    }
    
    ## for each bin on a chromosome
    #b = 1
    cat(paste0('Number of windows: ', length(Windows), '; window  sizes: ', paste(sapply(Windows, length), collapse=':'), '\n'))
    
    #if(any(sapply(Windows, length)<300)) stop('Number of SNPs in bins <300, nSNPs should be lowered')
    #b <- 1
    out[[i]] <- mclapply(1:length(Windows), function(b){
      cat(paste('Working on chromosome', i ,', window', b, '\n'))
      Window <- Windows[[b]]
      snp.id.bin <- snp.id[which(snp.id%in%Window)]
      snp_bin  <- snp[,which(snp.id%in%Window)]
      
      
      snpgdsCreateGeno(paste0(b,"_data2.gds"), genmat = snp_bin, sample.id = 1: nrow(snp) , snp.id = snp.id.bin, snpfirstdim=FALSE) ### sample ID should be 1: how many individuals you have
      file <- snpgdsOpen(paste0(b,"_data2.gds"))
      
      # estimate LD within a sliding window
      if(w2>length(snp.id.bin)) w2 = -1
      
      MAT <- snpgdsLDMat(file, method="r",slide=w2, verbose=FALSE)
      snpgdsClose(file)
      system(paste0('rm ', paste0(b,"_data2.gds")))
      
      LDmat <- data.table(remove=rep(NA, length(MAT$snp.id)))
      for(x in 1:(length(MAT$snp.id))){
        LDmat[(1:w2+x)[!is.nan(MAT$LD[,x])],as.character(x) := MAT$LD[!is.nan(MAT$LD[,x]),x]^2]
      }
      
      LDmat <- as.matrix(LDmat)
      LDmat <- LDmat[,-1]
      LDmat[is.na(LDmat)]<- 0
      LDmat[upper.tri(LDmat)] <- NA
      diag(LDmat) <- NA
      colnames(LDmat) <- snp.id.bin
      rownames(LDmat) <- snp.id.bin
      
      
      #cat(LDmat[1:10, 1:10])
      
      ## get raw data for LDna, using 'complete' linkage clustering and extracting mean and minium LD for each cluster (and its parent)
      if(verbose) cat('Preparing data for LDna \n')
      
      ldna <- LDnaRaw(LDmat, digits = 5, method='complete', fun=function(x){c(median(x, na.rm=TRUE), min(x, na.rm=TRUE))})
      
      if(verbose) cat('Extracting clusters \n')
      
      ## extract clusters based on teh extract.fun and the two LD thresholds
      out <- extractClusters(ldna, 
                             min.edges=1, 
                             extract = TRUE, 
                             rm.COCs =F, 
                             LD_threshold1 = LD_threshold1, 
                             LD_threshold2=LD_threshold2, 
                             plot.tree=plot.tree, 
                             plot.graph=FALSE, 
                             extract.fun=function(ldna, LD_threshold1, LD_threshold2){
                               colnames(ldna$clusterfile)[
                                 c(0,ldna$lambda_min$V2)>=LD_threshold1 & c(0,ldna$lambda_min$V4)<LD_threshold1 &
                                   c(0,ldna$lambda_min$V1)>=LD_threshold2] 
                             })
      
      clusters <- out[[1]]
      nested <- out[[2]]
      
      dimnames(nested) <- list(names(clusters), names(clusters))
      
      temp <- ldna$clusterfile[,colnames(ldna$clusterfile) %in% names(clusters)]
      
      ## remove any SOCs if its parent SOC also meets threshold requirements
      clusters <- clusters[!names(clusters) %in% colnames(temp)[which(nested=='COC', arr.ind = T)[,2]]]
      
      stats <- do.call(rbind, lapply(clusters, function(x){
        loci <- as.vector(unlist(x))
        temp2 <- LDmat[which(rownames(LDmat) %in% loci), which(colnames(LDmat) %in% loci)]
        c(Median=signif(median(na.omit(as.vector(temp2))), digits=5), MAD=signif(mad(as.vector(temp2), na.rm=TRUE, constant=1),3))
      }))
      
      
      ####
      if(any(!snp.id.bin %in% unlist(clusters))){
        ldna2 <- LDnaRaw(LDmat[!snp.id.bin %in% unlist(clusters), !snp.id.bin %in% unlist(clusters)], digits = 5, method='complete', fun=function(x){c(median(x, na.rm=TRUE), min(x, na.rm=TRUE))})
        
        ## extract clusters based on teh extract.fun and the two LD thresholds
        out <- extractClusters(ldna2, 
                               min.edges=1, 
                               extract = TRUE, 
                               rm.COCs =F, 
                               LD_threshold1 = LD_threshold2, 
                               LD_threshold2= LD_threshold2, 
                               plot.tree=plot.tree, 
                               plot.graph=FALSE, 
                               extract.fun=function(ldna2, LD_threshold1, LD_threshold2){
                                 colnames(ldna2$clusterfile)[
                                   c(0,ldna2$lambda_min$V2)>=LD_threshold2 & c(0,ldna2$lambda_min$V4)<LD_threshold2
                                   ] 
                               })
        
        clusters2 <- out[[1]]
        nested2 <- out[[2]]
        
        if(!is.null(clusters2)){
          
          if(!is.null(nested2)){
            dimnames(nested2) <- list(names(clusters2), names(clusters2))
            
            temp <- ldna2$clusterfile[,colnames(ldna2$clusterfile) %in% names(clusters2)]
            
            ## remove any SOCs if its parent SOC also meets threshold requirements
            clusters2 <- clusters2[!names(clusters2) %in% colnames(temp)[which(nested2=='COC', arr.ind = T)[,2]]]
          }
          
          clusters <- c(clusters, clusters2)
        }
        
        
        
        stats <- rbind(stats, do.call(rbind, lapply(clusters2, function(x){
          loci <- as.vector(unlist(x))
          temp2 <- LDmat[which(rownames(LDmat) %in% loci), which(colnames(LDmat) %in% loci)]
          c(Median=signif(median(na.omit(as.vector(temp2))), digits=5), MAD=signif(mad(as.vector(temp2), na.rm=TRUE, constant=1),3))
        })))
        ## merge clusters with non-clustered loci
        
        stats <- rbind(stats, cbind(rep(1, length(which(!snp.id.bin %in% unlist(clusters)))),rep(0, length(which(!snp.id.bin %in% unlist(clusters))))))
        
        clusters <- c(clusters, as.list(snp.id.bin[!snp.id.bin %in% unlist(clusters)]))
        
      }
      
      if(!is.null(plot.network)){
        ## plot network 
        cat('plotting network \n')
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
        
        nClust <- length(which(sapply(clusters, length)>1))
        
        temp1 <- sample(rep(col_vector, ceiling(nClust/length((col_vector)))))
        
        Col <- unlist(lapply(1:length(clusters), function(x) {
          ifelse(sapply(clusters, length)[x]==1, return('black'),  return(rep(temp1[1:nClust][x], sapply(clusters, length)[x])))
        }))[match(colnames(LDmat), unlist(clusters) )]
        
        
        temp2 <- sample(rep(col_vector[-1], ceiling(nClust/length((col_vector)))))
        
        frame.col <- unlist(lapply(1:length(clusters), function(x) {
          ifelse(sapply(clusters, length)[x]==1, return('black'),  return(rep(temp2[1:nClust][x], sapply(clusters, length)[x])))
        }))[match(colnames(LDmat), unlist(clusters) )]
        
        #plot.netork <- 'network'
        #threshold_net <- 0.9
        png(paste0(paste(plot.network, i, b, sep='_'), '.png'), res = 300, width = 6, height=6, units = 'in')
        plotLDnetwork(ldna,LDmat,option=1, threshold=threshold_net, col=Col, frame.color = frame.col)
        dev.off()
      }
      
      
      
      ## extract as many PCs from each cluster as necessary to explain 80% of the variation
      
      if(verbose) cat('Extracting PCs \n')
      
      PCs <- lapply(clusters, function(cl){
        ifelse(length(cl)>1, return(PC_score(snp_bin[,snp.id.bin%in%cl], PC_threshold)), return(snp_bin[,snp.id.bin==cl]))
      })
      
      
      temp <- lapply(PCs, function(x) ifelse(is.matrix(x), return(colnames(x)), return(NULL)))
      
      cat(paste0('PCs explaining ', round(100*mean(as.numeric(sapply(strsplit(sapply(temp[!sapply(temp, is.null)], function(x) ifelse(length(x)==2, return(x[2]), return(x[1]))), ':'), function(k) k[2])))), '% of the variation \n' ))
      
      temp2 <- table(sapply(PCs, function(x) ifelse(is.matrix(x), ncol(x), 1)))
      
      if(length(temp2)==1) temp2 <- c(temp2, 0); names(temp2) <- c('1', '2')
      
      cat(paste0('Reduction of number of tests: ', signif((1-sum(as.numeric(names(temp2))*temp2)/length(unlist(clusters)))*100, 4), '% \n'))
      
      return(list(reduction=
                    1-sum(as.numeric(names(temp2))*temp2)/length(unlist(clusters))
                  , clusters=clusters, PCs=PCs, stats=stats))
    }, mc.cores=mc.cores)
    
    #, mc.cores=mc.cores)}
  }
  
  ### prepare output
  cat('preparing output \n')
  if(length(out)==1)out <- c(list(NULL), out)
  out <- out[sapply(out, length)>0]
  
  
  temp <- lapply(out, function(chromosome){
    lapply(chromosome, function(window){
      window[[2]]
    })
  })
  
  bin_snps <- lapply(flatten(temp), as.numeric)
  
  # PCs
  temp <- flatten(lapply(out, function(chromosome){
    lapply(chromosome, function(window){
      window[[3]]
    })
  }))
  
  
  PVE <- signif(as.numeric(unlist(lapply(temp, function(x) {
    ifelse(is.matrix(x), return(sapply(strsplit(colnames(x), ':'), function(y) y[2])), return(1))
  }))), 3)
  
  
  bin_PCs <- lapply(temp, as.matrix)
  
  ## create map data ##
  
  Chr <- unlist(lapply(1:length(out), function(chromosome){
    lapply(1:length(out[[chromosome]]), function(window){
      rep(unique(mapsnp[,1])[chromosome], ncol(do.call(cbind, out[[chromosome]][[window]][[3]])))
    })
  }))
  
  Window <- unlist(lapply(1:length(out), function(chromosome){
    lapply(1:length(out[[chromosome]]), function(window){
      rep(window, ncol(do.call(cbind, out[[chromosome]][[window]][[3]])))
    })
  }))
  
  # chromosome <- out[[1]]
  # window <- chromosome[[1]]
  # x <- window[[2]][[1]]
  Pos <-  unlist(sapply(out, function(chromosome){
    sapply(chromosome, function(window){
      rep(sapply(window[[2]], function(x) mean(mapsnp[as.numeric(x),2])), sapply(window[[3]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)))
    })
  }))
  
  PC <- unlist(sapply(out, function(chromosome){
    sapply(chromosome, function(window){
      sapply(sapply(window[[3]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)), function(bin) 1:bin)
    })
  }))
  
  nSNPs <-  unlist(sapply(out, function(chromosome){
    sapply(chromosome, function(window){
      rep(sapply(window[[2]], length), sapply(window[[3]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)))
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
      rep(window[[4]][,1], sapply(window[[3]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)))
    })
  }))
  
  MAD <-  unlist(sapply(out, function(chromosome){
    sapply(chromosome, function(window){
      rep(window[[4]][,2], sapply(window[[3]], function(bin) ifelse(is.matrix(bin), ncol(bin), 1)))
    })
  }))
  
  Min <-  sapply(bin_snps, function(x) min(mapsnp[x,2]))[Grp]
  Max <-  sapply(bin_snps, function(x) max(mapsnp[x,2]))[Grp]
  Range <- abs(Min-Max)
  
  
  mapbin <- data.frame(data.table(Chr, Window, Pos=Pos, Min, Max, Range, nSNPs=nSNPs, Median=Median, MAD=MAD, PC=PC, PVE, Grp))
  colnames(mapbin) <- c('Chr', 'Window', 'Pos', 'Min', 'Max', 'Range', 'nSNPs', 'Median', 'MAD', 'PC', 'PVE', 'Grp')
  
  PCbin <- do.call(cbind, bin_PCs)
  
  ## order data by position
  # PCbin=PCbin[,order(mapbin$Chr, mapbin$Window, mapbin$Pos)]
  # 
  # temp <- mapbin[PC==1,]
  # 
  # bin_snps=bin_snps[order(temp$Chr,temp$Window, temp$Pos)]
  # 
  # mapbin=mapbin[order(mapbin$Chr, mapbin$Window, mapbin$Pos),] 
  # 
  return(list(mapbin=mapbin, PCbin=PCbin, bin_snps=bin_snps))
  
}

PC_score <- function(x_A, PC_threshold){
  cMean_A <- colMeans(x_A)
  for (i in 1:dim(x_A)[2]) x_A[, i] <- x_A[, i] - cMean_A[i]
  eigen_A <- eigen(var(x_A))
  scores_A <- x_A %*% eigen_A$vectors
  percentage_A = (cumsum(eigen_A$values))/sum(eigen_A$values)
  a = percentage_A >= PC_threshold
  a[which(a)[-1]] <-FALSE
  a[1] = TRUE
  zeta = as.matrix(scores_A[, a])
  colnames(zeta) <- paste(1:length(percentage_A[a]), percentage_A[a], sep=':')
  return(zeta)
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
