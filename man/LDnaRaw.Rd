\name{LDnaRaw}
\alias{LDnaRaw}
\title{Produces necessary data for LDna from a matrix of pairwise LD values}
\usage{
LDnaRaw(LDmat)
}
\arguments{
  \item{LDmat}{Lower diagonal matrix of pairwise LD values}
}
\value{
A list with two objects: \code{clusterfile} and
\code{stats}. The first, \code{clusterfile} is a matrix
with all unique clusters (defined by a specific set of
loci) as columns and each locus as rows. In this matrix
\code{TRUE} indicates the presense of a locus in a specific
cluster, else \code{FALSE} is given. The second,
\code{stats} is a data frame that gives all the edges for
the \code{single linkage clustering tree} along with the
following information for each edge: probe, parent probe,
and the corresponding \code{nV}, \code{nE}, and the lambda
for each probe.
}
\description{
Takes a matrix of pairwise LD values and produces data
files for subsequent LD network analyses.
}
\details{
In LDna, trees are used to visualise the process of cluster
merging with decreasing LD threshold, where branches
represent clusters/loci and branch points represent the
merger of one or several clusters/loci. Distance indicates
LD threshold values at which these events occur. To
identify mergers of distinct clusters, the change in median
LD between all loci in a cluster, before and after merger,
lambda, is recorded for all clusters in the tree.  With a
lower diagonal matrix of LD values, the LD values are
transformed into dissimilarities (1-LD) and a \code{single
linkage clustering} algorithim (see \code{\link{hclust}})
is used to produce a tree that describes the order of
cluster merging with decreasing LD threshold. For each
cluster the number of vertices, \code{nV}, the number of
edges, \code{nE} are recorded.
}
\examples{
# Simple upper diagonal LD matrix
LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.44, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.36, 0.24, NA, NA, NA, 0.48, 0.32, 0.2, 0.36, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.24, NA, NA, NA, NA, NA, NA, 0.44, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
# Calculate raw data
ldna <- LDnaRaw(LDmat)
# For larger data set
data(LDna)
ldna <- LDnaRaw(r2.baimaii_subs)
}
\author{
Petri Kemppainen \email{petrikemppainen2@gmail.com},
Christopher Knight \email{Chris.Knight@manchester.ac.uk}
}
\keyword{LDnaRaw}

