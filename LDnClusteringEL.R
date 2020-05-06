
LDnClusteringEL <- function(EL_path = "./LD_EL", nSNPs=1000, w1=10, w2=100, LD_column = 4, out_folder="LDnCl_out", 
                               LD_threshold1 = 0.5, LD_threshold2 = 0.8, PC_threshold=0.8, mc.cores=1, plot.network=NULL, threshold_net=0.9, to_do = NULL, min.cl.size=2){
  
  out <- list()
  EL_files <- list.files(EL_path)[to_do]
  temp <- list.files()
  if(!any(grepl(out_folder, temp))) system(paste("mkdir ",out_folder))
  
  # i <- 1
  for(i in 1:length(EL_files)){
    cat(paste0("Finding windows in file: ", EL_files[i], "\n"))
    
    chr <- gsub(pattern = "group", replacement = "", x = strsplit(EL_files[i], "_")[[1]][2])
    
    el <- fread(paste(EL_path, EL_files[i], sep = "/"))
    ## get ld between adjacent pairs of loci
    setnames(el, colnames(el)[c(1,2, LD_column)], c("from", "to", "r2"))
    
    el[,from := as.character(from)]
    el[,to := as.character(to)]
    #Columns <- c(1,2,LD_column)
    #dim(el)
    el_adj <- as.matrix(el[!duplicated(from),.(from, to, r2)])
    snp.pos <- as.vector(unlist(el_adj[,1]))
    snp.pos <- c(snp.pos, as.vector(unlist(el_adj[nrow(el_adj),2])))
    unique(c(el_adj[,1], el_adj[,2])) == snp.pos
    #dim(el_adj)
    #snp.pos <- as.vector(unique(unlist(el_adj[,1:2])))
    
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
    #308*200*20
    # w <- 2
    # i <- 1
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
      #V(g)$name
      d_g <- decompose.graph(g)
      d_g <- d_g[sapply(d_g, function(x) length(V(x)$name)>1)]
      
      #sapply(d_g, function(x) length(V(x)$name))
      
      #snp_w <- matrix(sample(c(0,1,2), length(Window)*100,replace = TRUE), 100, length(snp.pos_w))
      
      LDmat <- get.adjacency(g,attr = 'weight',names = TRUE,sparse = FALSE)
      diag(LDmat) <- NA
      LDmat <- t(LDmat)
      LDmat[upper.tri(LDmat)] <- NA
      #LDmat[1:20, 1:20]
      #clade <- tree$edge[,2][tree$edge[,1]==tree$edge[1,1]]
      #x <- clade[1]
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
            #& length(cl)>2 | Median.temp > LD_threshold1 & length(cl)==2
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
            # }else{
            #   clusters.sub[[x]] <<- paste(chr, cl, sep=":")
            #   MCL.sub[[x]]  <<-cl[which.max(apply(LDmat.part, 1, function(h) median(h, na.rm = TRUE)))]
            #   Median.sub[[x]] <<- median(LDmat.part, na.rm  = TRUE)
            #   Mad.sub[[x]] <<- mad(LDmat.part, na.rm  = TRUE)
            # }
          }
        })
      }
      # snp_w[,which]
      # dim(snp_w)
      # length(snp.pos_w)
      # clusters.sub[[x]] <- paste(chr, cl, sep=":")
      # PC_sc <- PC_score(snp_w[,snp.pos_w %in% cl], PC_threshold)
      # PVE.sub[[x]] <<- PC_sc[[2]]
      # PCs.sub[[x]] <<- PC_sc[[1]]
      # MCL.sub[[x]] <<- snp_cl[,which.max(apply(LDmat.part, 1, function(h) median(h, na.rm = TRUE)))]
      # Median.sub[[x]] <<- Median.temp
      # Mad.sub[[x]] <<- mad(LDmat.part, na.rm  = TRUE)
      # 
      
      clusters <- list()
      MCL <- list()
      Median <- list()
      Mad <- list()
      
      for(y in 1:length(d_g)){
        #print(y)
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
      
      
      # clusters <- flatten(clusters)
      # PVE <- flatten(PVE)
      # PCs <- flatten(PCs)
      # MCL <- flatten(MCL)
      # Cons <- flatten(Cons)
      # Mad <- unlist(flatten(Mad))
      # Median <- unlist(flatten(Median))
      # 
      clusters <- flatten(clusters)
      #sum(table(sapply(clusters, function(x) length(x)))[-1])
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
    }, mc.cores=4)
    #length(out)
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
