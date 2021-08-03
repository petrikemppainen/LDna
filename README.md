R-package: LDna
-------------
Linkage disequilibrium (LD) network analysis (**LDna**) can be used to find clusters of loci in high LD from population genomic data sets using network analysis. It provides a means to partition population genomics data into sets of loci that bear similar population genetic signals. It can for instance be used to *in silico* identify inversion polymorphism and loci involved in local adaptation. As it only requires a matrix of pair-wise LD values it is particularly useful for non-model species where closely related and well-characterized reference genomes are not available.

LDna also includes a function that can be used to find cluster of loci connected by high LD along non-overlapping windows and uses Principal Component Analysis to summerise the information from the clusters, to increase power (less conservative multiple correction) in Genome Wide Association Studies (GWAS).

Reference: http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12369/abstract

For this we recommend the *LDnClusteringEL* which requires edge-list of pairwise LD-values as input. For next-gen data we recommend using **ngsLD** (https://github.com/fgvieira/ngsLD) for LD-estimation.

LDna v2.0 contains some major upgrades (see full version history below) mainly the inclusion of two new functions (*LDnaRaw2* and *extractBranches*) that allow for much simpler detection of LD-clusters. The old versions (*LDnaRaw* and *extractClusters*) are still included as the tutorials for these are much more comprehensive and they still work well for smaller data sets and allow for much better visualisation of the principles of **LDna**. For this you might however be better off by using the v0.64 instead.

The function *extractBranches* is a much simplified version of *extractClusters* that only relies on the parameter *|E|min* (the minimum number of edges for LD-clusters) which essentially determines how many branches are allowed in a tree. Then all branches are considered as LD-clusters, thus lower *|E|min* values lead to many smaller and more strongly correlated LD-clusters and conversely higher *|E|min* values lead to fewer but larger LD-clusters where loci will on average be less correlated. Smaller clusters potentially have stronger correlations with traits of interest (in GWAS) but also more conservative corrections for multiple testing and vice versa.

For more details please see https://www.biorxiv.org/content/10.1101/2021.01.26.428263v1, in particular Supporting Methods 4 to understand the trade off between high/low parameter settings for *|E|min*. Please note that in this paper LDna-1 is equivalent to *LDnClusteringEL*

For questions, email me at petrikemppainen2@gmail.com

For small data sets an getting famliar with LDna basics, I still recommend checking out v.0.64!

###Installing

With **devtools** (accessible from CRAN) **LDna** can be installed by:
```r
devtools::install_github("petrikemppainen/LDna", ref = 'v2.0')
```
This downloads the source directly from **github** and builds the vignettes and thus requires LaTeX to be installed on your computer.

Alternatively, download the source file (LDna_v2.0.tar.gz) directly and install by:
```r
install.packages("/path_to/source_file", repos = NULL, type = "source")
```
Please install and follow package documentation (includes two tutorials) for more information.

###Version history

v.63 was used for LDnClustering (https://www.ncbi.nlm.nih.gov/pubmed/29673105).

There's an implementation that takes a path to a folder with vcf files (one for each chromosome) as input (LDnClusteringVCF). Both this and LDnClustering (which requires a snp and a map file) now also choose one SNP from each cluster based on either the locus with highest intra cluster median LD or creates a consensus SNP based on all the SNPs from the cluster. These can be anlysed with any GWAS methods that take SNP data as input. As long as LD thresholds are high (i.e. SNPs are strongly correlated within each cluster) this should work fine. With low threshold settings for LD the PC method is expected to better capture the information from each cluster, but this has not been tested properly yet.

There's some example code on how to use emmax for LD-clustered data. Tutorials are still in progress. 

LDnClustering now removes all edges below LD_threshold1 and decomposes this graph (essentially performing single linkage clustering) and in the resulting sub-clusters LDnClustering performs complete linkage clustering and recursively finds clades where median LD is below LD_threshold2. In this way the hireachichal clustering does not have to be done on the whole tree (window) but in the sub-clusters instead. Also, the recursive step to find relevant clades means that less singleton clusters are produced compared to the original implementation.

As of v.653 I have added the possibility multithreading for *LDnClusteringVCF* (where a folder with "vcf.gz" files can be specified, one file for each chromosome/linkage group) and improved speed. We also added the possibility to use pre-calculated pairwise LD values in the form of edge lists (first column=locus 1, second column=locus 2, and the column where the  r^2 value for this pairwise comparison can be given by *LD_column*) as input (*LDnClusteringEL*). We use this for instance when we estimate r^2 from genotype likelihoods.

As of v.7 I have added "branch traversal" to function "extractClusters", Branch traversal means that as long as a branch contains outliers, the outlier clusters (SOC or COC) that is closest to the base of the branch is selected as the outlier cluster. This ensures that a wide range threshold values for \eqn{\lambda} gives similar results. 

This version does no longer contain the original vignettes as they don't work anymore with the updated code. See examples for each function instead. The vignettes are still aveailable from prevous versions.

