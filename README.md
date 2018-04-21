R-package: LDna
-------------

### Note: While v.63 was used for LDnClustering and will remain unchanged (also available from Dryad) and works well LDnClustering, I had to temporarily compromise some of the original functionality of LDna. Here I will attempt to fix this, as well as update documentation and produce some example code and tutorials for LDnClustering. I will also try to implement a much faster algorithim such that LDnClustering can more readily be used for whole genome data sets.

### Updatae, the faster implementation is now functional, and there's some example code on how to use emmax for LD-clustered data for the function 'LDnClustering'. Tutorials are still in progress. Rather than performing two steps of comlete linkage clustering, LDnClustering now removes all edges below LD_threshold1 and decomposes this graph (essentially performing single linkage clusetring) and in the resulting sub-clusters LDnClustering performs complte linkage clustering and recursively finds clades where median LD is below LD_threshold2. In this way the hireachichal clustering does not have to be done on the whole tree but in the sub-clusters instead. Also the recursive step to find relevant clades means that less singleton clusters are produced compared to the original implementation.


Linkage disequilibrium (LD) network analysis (LDna) can be used to find clusters of loci in high LD from population genomic data sets using network analysis. It provides a means to partition population genomics data into sets of loci that bear similar population genetic signals. It can for instance be used to *in silico* identify inversion polymorphism and loci involved in local adaptation. As it only requires a matrix of pair-wise LD values it is particularly useful for non-model species where closely related and well-characterized reference genomes are not available.

LDna also includes a function that can be used to find cluster of loci connected by high LD along non-overlapping windows and uses Principal Component Analysis to summerise the information from the clusters, to increase power (less conservative multiple correction) in Genome Wide Association Studies (GWAS). 
Reference: http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12369/abstract

Current beta version is 0.64.

Any questions or suggestions may be posted at: https://groups.google.com/forum/#!forum/ld-network-analysis

###Installing

With **devtools** (accessible from CRAN) **LDna** can be installed by:
```r
devtools::install_github("petrikemppainen/LDna", ref = 'v.64')
```
This downloads the source directly from **github** and builds the vignettes and thus requires LaTeX to be installed on your computer.

Alternatively, download the source file (LDna_0.64.tar.gz) directly and install by:
```r
install.packages("/path_to/source_file", repos = NULL, type = "source")
```
Please install and follow package documentation (includes two tutorials) for more information.
