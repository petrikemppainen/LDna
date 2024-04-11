R-package: LDna
-------------
Linkage disequilibrium (LD) network analysis (**LDna**) can be used to find clusters of loci in high LD from population genomic data sets using network analysis. It provides a means to partition population genomics data into sets of loci that bear similar population genetic signals. It can for instance be used to *in silico* identify inversion polymorphism and loci involved in local adaptation. As it only requires a matrix of pair-wise LD values it is particularly useful for non-model species where closely related and well-characterized reference genomes are not available.

LDna also includes a function that can be used to find cluster of loci connected by high LD along non-overlapping windows and uses Principal Component Analysis to summerise the information from the clusters, to increase power (less conservative multiple correction) in Genome Wide Association Studies (GWAS).

Reference: http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12369/abstract

For this we recommend the *LDnClusteringEL* which requires edge-list of pairwise LD-values as input. For next-gen data we recommend using **ngsLD** (https://github.com/fgvieira/ngsLD) for LD-estimation.

From LDna v2.0 *extractClusters* has been replaced with *extractBranches* which is a much simplified version of *extractClusters* that only relies on the parameter *|E|min* (the minimum number of edges for LD-clusters) which essentially determines how many branches are allowed in a tree. Then all branches are considered as LD-clusters, thus lower *|E|min* values lead to many smaller and more strongly correlated LD-clusters and conversely higher *|E|min* values lead to fewer but larger LD-clusters where loci will on average be less correlated. Smaller clusters potentially have stronger correlations with traits of interest (in GWAS) but also more conservative corrections for multiple testing and vice versa.

LDna v2.14 fixes some bugs in *extractBranches* that sometimes did not extract some of the LD-clusters that should have been extracted. Furthermore, in small data sets and using *min.edges=0* it also gives you all singleton clusters (i.e. those loci that are not part of any other LD-clusters) and includes them in the tree.

While this version does not have any manual or vignettes, pdf versions for v.63 can still be found at https://github.com/petrikemppainen/LDna/tree/v.63/inst/doc which will give you the basics. Please refer to the examples and information for each function for further details.

For more details please see https://www.biorxiv.org/content/10.1101/2021.01.26.428263v1, in particular Supporting Methods 4 to understand the trade off between high/low parameter settings for *|E|min*. Please note that in this paper LDna-1 is equivalent to *LDnClusteringEL*

Any questions or suggestions may be posted at: https://groups.google.com/forum/#!forum/ld-network-analysis

###Installing

With **devtools** (accessible from CRAN) **LDna** can be installed by:
```r
devtools::install_github("petrikemppainen/LDna", ref = 'v.2.14')
```
This downloads the source directly from **github** and builds the vignettes and thus requires LaTeX to be installed on your computer.

Alternatively, download the source file (LDna_v.2.14.tar.gz) directly and install by:
```r
install.packages("/path_to/source_file", repos = NULL, type = "source")
```
Please install and follow package documentation for more information.
