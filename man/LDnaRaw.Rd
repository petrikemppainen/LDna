\name{LDnaRaw}
\alias{LDnaRaw}
\title{Produces necessary data for LDna from a matrix of pairwise LD values}
\usage{
LDnaRaw(LDmat)
}
\arguments{
  \item{LDmat}{Lower diagonal matrix of pairwise LD values,
  \eqn{R^2} is strongly recommended.}
}
\value{
A list with two objects: \code{clusterfile} and
\code{stats}. \code{Clusterfile} is a matrix with all
unique clusters as columns and each locus as rows, where
\code{TRUE} indicates the presense of a locus in a specific
cluster (else \code{FALSE} is given). \code{Stats} is a
data frame that gives all edges for the \code{single
linkage clustering tree} along with the following
information for each edge: cluster (focal cluster), parent
cluster (cluster efter merger), and the corresponding
\code{nV}, \code{nE}, and \eqn{\lambda} for each cluster
(cluster).
}
\description{
Takes a lower diagonal matrix of pairwise LD values and
produces data files for subsequent LD network analyses.
}
\details{
In LDna clusters represent loci (vertices) connected by LD
values (edges) above given thresholds. From a lower
diagonal matrix of LD values a \code{single linkage
clustering algorithim} (see \code{\link{hclust}}) is used
to produce a tree that describes the order of cluster
merging with decreasing LD threshold. In these trees
branches represent clusters/loci and branch points
represent the merger of one or several clusters/loci and
distance indicates LD threshold values at which these
events occur. To identify mergers of distinct clusters, the
change in median LD before and after merger, \eqn{\lambda},
is calculated for all clusters in the tree. \eqn{\lambda}
is defined as \eqn{(x_1-x_2)*n} where \eqn{x_1} is the
median of all pairwise LD values within the focal cluster,
\eqn{x_2} is the median of all pairwise values between loci
within the focal cluster and between the loci in the focal
cluster and any other loci it merges with, and \eqn{n} is
the number of loci in the focal cluster. High values of
\eqn{\lambda}, relative to the median of all values across
the tree, indicate the merging of clusters with highly
distinct LD signals. For each cluster the number of
vertices, \code{nV}, the number of edges, \code{nE} are
recorded.
}
\examples{
# Simple upper diagonal LD matrix
LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.44, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.36, 0.24, NA, NA, NA, 0.48, 0.32, 0.2, 0.36, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.24, NA, NA, NA, NA, NA, NA, 0.44, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
# Produce raw data
ldna <- LDnaRaw(LDmat)
# Larger data set
data(LDna)
ldna <- LDnaRaw(r2.baimaii_subs)
}
\author{
Petri Kemppainen \email{petrikemppainen2@gmail.com},
Christopher Knight \email{Chris.Knight@manchester.ac.uk}
}
\seealso{
\code{\link{extractClusters}}, \code{\link{summaryLDna}}
and \code{\link{plotLDnetwork}}
}
\keyword{LDnaRaw}

