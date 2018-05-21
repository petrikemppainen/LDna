#' Combines data from files produced by \code{\link{LDnClusteringVCF}}
#'
#' \code{\link{LDnClusteringVCF}} produces one \code{.rds} file for each chromosome (from separate .vcf files). This function concatenates the information into a single file
#' 
#' @param LDnCl_out Path to folder containing output from \code{\link{LDnClusteringVCF}}
#' @keywords Linkage disequilibrium network clustering, complexity reduction
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}, zitong.li \email{lizitong1985@@gmail.com}
#' @return Returns concatenated file with the same data as in the individual \code{.rds} files. Can be used with e.g. \code{\link{emmax_group}} as with output produced by \code{\link{LDnClustering}}
#' @references Kemppainen, P., Knight, C. G., Sarma, D. K., Hlaing, T., Prakash, A., Maung Maung, Y. N., Walton, C. (2015). Linkage disequilibrium network analysis (LDna) gives a global view of chromosomal inversions, local adaptation and geographic structure. Molecular Ecology Resources, 15(5), 1031-1045. https://doi.org/10.1111/1755-0998.12369\cr
#' \cr
#' Li, Z., Kemppainen, P., Rastas, P., Merila, J. Linkage disequilibrium clustering-based approach for association mapping with tightly linked genome-wide data. Accepted to Molecular Ecology Resources.
#' @examples
#' \dontrun{
#' LDnCL_combined <- Combine_LDnCl(path/to/rds/files)
#'}                   
#' @export
#' 
Combine_LDnCl <- function(LDnCl_out){
  file_names <- list.files(LDnCl_out)
  
  LDnCl_combined <- list()
  LDnCl_combined$cluster_summary  <- do.call(rbind, lapply(file_names, function(x){
    readRDS(paste(LDnCl_out, x, sep = '/'))$cluster_summary
  }))
  
  LDnCl_combined$cluster_PCs <- do.call(cbind, lapply(file_names, function(x){
    readRDS(paste(LDnCl_out, x, sep = '/'))$cluster_PCs
  }))
  
  x <- file_names[1]
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
