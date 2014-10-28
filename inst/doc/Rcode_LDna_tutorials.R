### R code for LDna tutorials
## R code for basic tutorial ##############
library(LDna)
##
LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.44, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.36, 0.24, NA, NA, NA, 0.48, 0.32, 0.2, 0.36, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.24, NA, NA, NA, NA, NA, NA, 0.44, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(paste("L", 1:8, sep=""), paste("L", 1:8, sep="")))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
LDmat
##
ldna <- LDnaRaw(LDmat)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
ldna$clusterfile
ldna$stats
##
par(mfcol=c(1,2)) # this will print two figures side by side     
clusters1 <- extractClusters(ldna, min.edges = 0, phi = 1)
clusters1
clusters1[[1]]
clusters1[[2]]
##
clusters2 <- extractClusters(ldna, min.edges = 0, phi = 0.25, rm.COCs=FALSE)
clusters2
##
summary <- summaryLDna(ldna, clusters2, LDmat)
summary
##
plotLDnetwork(LDmat=LDmat, option=1, threshold=0.6)
plotLDnetwork(LDmat=LDmat, option=1, threshold=0.40)
###########
# Creates a 'weighted' graph object from our lower diagonal LD matrix
g <- graph.adjacency(LDmat, mode="lower", diag=FALSE, weighted=T)
# Removes edges with weights (LD values) below 0.5
g <- delete.edges(g, which(E(g)$weight<=0.5))
# Removes 'unconnected' vertices (loci)
g <- delete.vertices(g, which(degree(g) < 1))
# plots graph, 'edge.label' is used to add the LD values for each edge
plot.igraph(g, edge.label=E(g)$weight)
###########
par(mfcol=c(1,3))
plotLDnetwork(ldna, LDmat, option=2, clusters=clusters2, summary=summary)
plotLDnetwork(ldna, LDmat, option=2, clusters=clusters2, summary=summary, after.merger=TRUE)
#
clusters3 <- extractClusters(ldna, min.edges = 0, lambda.lim=0.5, rm.COCs=FALSE)
## R code for advnaced tutorial ###################
data(LDna)
LDmat <- r2.baimaii_subs
##
par(mfcol=c(1,3))
plotLDnetwork(LDmat=LDmat, option=1, threshold=0.8)
plotLDnetwork(LDmat=LDmat, option=1, threshold=0.5)
plotLDnetwork(LDmat=LDmat, option=1, threshold=0.3)
##
ldna <- LDnaRaw(LDmat)
##
par(mfcol=c(1,3))
extractClusters(ldna, min.edges = 1, plot.tree = TRUE, extract=FALSE)
extractClusters(ldna, min.edges = 5, plot.tree = TRUE, extract=FALSE)
extractClusters(ldna, min.edges = 18, plot.tree = TRUE, extract=FALSE)
##
par(mfcol=c(1,2))
plotLDnetwork(LDmat=LDmat, option=1, threshold=0.31)
plotLDnetwork(LDmat=LDmat, option=1, threshold=0.30)
##
par(mfcol=c(1,2))
clusters1 <- extractClusters(ldna, min.edges = 18, phi = 4, rm.COCs = FALSE, plot.tree = TRUE, plot.graph = FALSE)
summary1 <- summaryLDna(ldna, clusters1, LDmat)
summary1
##
par(mfcol=c(2,3))
plotLDnetwork(ldna, LDmat, option=2, clusters=clusters1, summary=summary1)
##
par(mfcol=c(1,2))
clusters_phi2 <- extractClusters(ldna, min.edges = 18, phi = 2, rm.COCs = FALSE, plot.tree = TRUE, plot.graph = TRUE)
summary_phi2 <- summaryLDna(ldna, clusters_phi2[names(clusters_phi2) %in% c("90_0.47","148_0.32")], LDmat)
summary_phi2
##
par(mfcol=c(1,2))
plotLDnetwork(ldna, LDmat, option=2, clusters=clusters_phi2[names(clusters_phi2) %in% c("90_0.47","148_0.32")], summary=summary_phi2, full.network=FALSE, include.parent=TRUE, after.merger=FALSE)
##
par(mfcol=c(1,2))
plotLDnetwork(ldna, LDmat, option=2, clusters=clusters_phi2[names(clusters_phi2) %in% c("90_0.47","148_0.32")], summary=summary_phi2, full.network=FALSE, include.parent=TRUE, after.merger=TRUE)
clusters2 <- extractClusters(ldna, min.edges=8, phi=2, rm.COCs=TRUE, plot.tree=TRUE, plot.graph=FALSE)
clusters2$"194_0.21"
##
par(mfcol=c(1,2))
summary2 <- summaryLDna(ldna, clusters=clusters2[names(clusters2) %in% c("118_0.38","194_0.21")], LDmat)
plotLDnetwork(ldna, LDmat, option=2, clusters=clusters2[names(clusters2) %in% c("118_0.38","194_0.21")], summary=summary2)
##
par(mfcol=c(1,2))
clusters_final <- extractClusters(ldna, min.edges = 18, phi=4)
clusters_final
summary_final <- summaryLDna(ldna, clusters_final, LDmat)
summary_final
##
par(mfcol=c(2,2))
plotLDnetwork(ldna, LDmat, option=2, clusters=clusters_final, summary=summary_final)
############### extra tutorial ################
# here we use file "r2.sim". It should already exist if you executed "load(LDna)" above
# a few functions to explore this file
is(r2.sim)
length(r2.sim)
r2.sim[[1]][1:10, 1:10]
is(r2.sim[[1]])
names(r2.sim)

# get ldna raw data for all data sets in the "r2.sim"
ldna <- list()
for(i in 1:length(r2.sim)){
  ldna[[i]] <- LDnaRaw(r2.sim[[i]])
}

i <- 2 # choose which data set to look at
extractClusters(ldna[[i]], extract=F)   # we are mainly intereste in the slink, so use 'extract=F'
plotLDnetwork(LDmat=r2.sim[[i]], option=1, threshold = 0.75); title(xlab=names(r2.sim)[[i]]) # change the threshold to see what happens

