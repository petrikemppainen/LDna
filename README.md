R-package: LDna
-------------

### Note: this is a beta version 0.63 which was used for 'Linkage disequilibrium cluster-based approach for association mapping with tightly linked genome-wide data' in review for Mol Ecol Res. 

It works well with code provided with this publication (available from dryad), but has not properly been tested with other data. Documentation is also not yet complete, I aim to add a tutorial specifically for the LDn-binning approach, in the near future. I also need to make sure existing users of LDna can continue using it as before.

Linkage disequilibrium (LD) network analysis (LDna) can be used to find clusters of loci in high LD from population genomic data sets using network analysis. It provides a means to partition population genomics data into sets of loci that bear similar population genetic signals. It can for instance be used to *in silico* identify inversion polymorphism and loci involved in local adaptation. As it only requires a matrix of pair-wise LD values it is particularly useful for non-model species where closely related and well-characterized reference genomes are not available.

Reference: http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12369/abstract

Current beta version is 0.63.

Any questions or suggestions may be posted at: https://groups.google.com/forum/#!forum/ld-network-analysis

###Installing

With **devtools** (accessible from CRAN) **LDna** can be installed by:
```r
devtools::install_github("petrikemppainen/LDna", ref = 'v.63')
```
This downloads the source directly from **github** and builds the vignettes and thus requires LaTeX to be installed on your computer.

Alternatively, download the source file (LDna_0.63.tar.gz) directly and install by:
```r
install.packages("/path_to/source_file", repos = NULL, type = "source")
```
Please install and follow package documentation (includes two tutorials) for more information.
