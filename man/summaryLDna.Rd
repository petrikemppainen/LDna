\name{summaryLDna}
\alias{summaryLDna}
\title{Summerizes LD clusters}
\usage{
summaryLDna(ldna, clusters, LDmat)
}
\arguments{
  \item{ldna}{Output from \code{\link{LDnaRaw}}}

  \item{clusters}{A file created by
  \code{\link{extractClusters}} that contains a list of all
  extracted clusters and their loci}

  \item{LDmat}{Lower diagonal matrix of pairwise LD values}
}
\value{
Returns a data frame with each row corresponding to a
cluster, in decreasing order with respect to highest LD
threshold value at which they are present. Column
\code{Type} specifies if a cluster is a "single outlier
cluster" (SOC) or a "compound oulier cluster" (COC). Column
"Merge.at" specfies the LD threshold at which the cluster
merges. Column "Mean.LD" gives the average LD of all
pairwise LD values between all loci in a cluster. The
remaining columns are self exlanatory.
}
\description{
Creates a data frame with summary information of LD
clusters extracted by the function
\code{\link{extractClusters}}
}
\examples{
# Simple upper diagonal LD matrix
LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.12, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.08, 0.4, NA, NA, NA, 0.48, 0.32, 0.04, 0.44, 0.36, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.32, NA, NA, NA, NA, NA, 0.72, 0.68, 0.28, NA, NA, NA, NA, NA, NA, 0.2, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
# Transform
# Get necessary data from  LD matrix
ldna <- LDnaRaw(LDmat)
# extract clusters and plot graphs, for this small example min.edges is set to zero such that each tip clade corresponds to an individual locus
clusters <- extractClusters(ldna, min.edges=0, phi=1, rm.COCs=FALSE)
summary <- summaryLDna(ldna, clusters, LDmat)
}
\author{
Petri Kemppainen \email{petrikemppainen2@gmail.com}
}
\seealso{
\code{\link{extractClusters}}, \code{\link{LDnaRaw}}
}
\keyword{summaryLDna}

