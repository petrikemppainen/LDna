\name{summaryLDna}
\alias{summaryLDna}
\title{Summerises LD clusters}
\usage{
summaryLDna(ldna, clusters, LDmat)
}
\arguments{
  \item{ldna}{Output from \code{\link{LDnaRaw}}}

  \item{clusters}{Output from
  \code{\link{extractClusters}}}

  \item{LDmat}{Lower diagonal matrix of pairwise LD values}
}
\value{
Returns a data frame with each row corresponding to a
cluster, in decreasing order with respect to highest LD
threshold value at which they are present. Column
\emph{Type} specifies if a cluster is a \emph{"single
outlier cluster", SOC} or a \emph{"compound oulier
cluster", COC}. Column \emph{"Merge.at"} specfies the LD
threshold for cluster merger. "Median.LD" gives the average
LD of all pairwise values between loci in a cluster, and
"MAD.LD" gives the scaled median absolute deviation for
these LD values.
}
\description{
Creates a data frame with summary information of LD
clusters extracted by function
\code{\link{extractClusters}}
}
\examples{
# Simple upper diagonal LD matrix
LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.44, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.36, 0.24, NA, NA, NA, 0.48, 0.32, 0.2, 0.36, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.24, NA, NA, NA, NA, NA, NA, 0.44, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
# Produce neccesary data for LDna
ldna <- LDnaRaw(LDmat)
# Extract clusters and plot graphs, for this small example min.edges is set to zero such that each tip clade corresponds to an individual locus
clusters <- extractClusters(ldna, min.edges=0, phi=0.25, rm.COCs=FALSE)
summary <- summaryLDna(ldna, clusters, LDmat)
}
\author{
Petri Kemppainen \email{petrikemppainen2@gmail.com}
}
\seealso{
\code{\link{extractClusters}}, \code{\link{LDnaRaw}} and
\code{\link{plotLDnetwork}}
}
\keyword{summaryLDna}

