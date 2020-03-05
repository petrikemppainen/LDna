R-package: LDna
-------------
Linkage disequilibrium (LD) network analysis (LDna) can be used to find clusters of loci in high LD from population genomic data sets using network analysis. It provides a means to partition population genomics data into sets of loci that bear similar population genetic signals. It can for instance be used to *in silico* identify inversion polymorphism and loci involved in local adaptation. As it only requires a matrix of pair-wise LD values it is particularly useful for non-model species where closely related and well-characterized reference genomes are not available.

LDna also includes a function that can be used to find cluster of loci connected by high LD along non-overlapping windows and uses Principal Component Analysis to summerise the information from the clusters, to increase power (less conservative multiple correction) in Genome Wide Association Studies (GWAS). 
Reference: http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12369/abstract

v.63 was used for LDnClustering (https://www.ncbi.nlm.nih.gov/pubmed/29673105).

There's an implementaiton that takes a path to a folder with vcf files (one for each chromsome) as input (LDnClusteringVCF). Both this and LDnClustering (which requires a snp and a map file) now also choose one SNP from each cluster based on either the locus with highest intra cluster median LD or creates a consensus SNP based on all the SNPs from the cluster. These can be anlysed with any GWAS methods that take SNP data as imput. As long as LD thresholds are high (i.e. SNPs are strongly correlated within each cluster) this should work fine. With low threshold settings for LD the PC method is expected to better capture the information from each cluster, but this has not been tested properly yet.

There's some example code on how to use emmax for LD-clustered data. Tutorials are still in progress. 

LDnClustering now removes all edges below LD_threshold1 and decomposes this graph (essentially performing single linkage clusetring) and in the resulting sub-clusters LDnClustering performs complete linkage clustering and recursively finds clades where median LD is below LD_threshold2. In this way the hireachichal clustering does not have to be done on the whole tree (window) but in the sub-clusters instead. Also, the recursive step to find relevant clades means that less singleton clusters are produced compared to the original implementation.

As of v.0.65, there's also function (LDnClusteringEL) that performs LDn-clustering based on precalculated LD values (as edge lists). The path to folder where these edge lists (one for each chromosome) resides. There are examples for this.
-
Any questions or suggestions may be posted at: https://groups.google.com/forum/#!forum/ld-network-analysis

###Installing

With **devtools** (accessible from CRAN) **LDna** can be installed by:
```r
devtools::install_github("petrikemppainen/LDna", ref = 'master')
```
This downloads the source directly from **github** and builds the vignettes and thus requires LaTeX to be installed on your computer.

Alternatively, download the source file (LDna_0.65.tar.gz) directly and install by:
```r
install.packages("/path_to/source_file", repos = NULL, type = "source")
```
Please install and follow package documentation (includes two tutorials) for more information.
