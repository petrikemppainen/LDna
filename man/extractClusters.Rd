\name{extractClusters}
\alias{extractClusters}
\title{Extracts LD clusters.}
\usage{
extractClusters(ldna, min.edges = 0, phi = 1, lambda.lim = NULL,
  rm.COCs = TRUE, extract = TRUE, plot.tree = TRUE, plot.graph = TRUE)
}
\arguments{
  \item{ldna}{File of class \code{ldna} that contains raw
  data for LDna analyses}

  \item{min.edges}{Minimum number of edges for a cluster
  that is allowed in a tree}

  \item{phi}{Controls \eqn{\lambda_{lim}} which sets the
  threshold above which \eqn{\lambda} are considered as
  outliers. Default is two, values below this are not
  recommended.}

  \item{rm.COCs}{If TRUE (default), automatically removes
  "compound oulier clusters" (COCs).}

  \item{extract}{If TRUE (default), returns a list of
  cluster names. If "FALSE" only prints a tree with no
  "oulier clusters" (OCs) indicated.}

  \item{plot.tree}{If TRUE (default), plots tree. Has no
  effect if extract=FALSE.}

  \item{plot.graph}{If TRUE (default), plots all
  \eqn{\lambda} ordered from low to high and indicates
  which values are above \eqn{\lambda_{lim}}.}

  \item{lambda.lim}{If not NULL gives a fixed value for
  \eqn{\lambda_{lim}}. Overrides \code{phi}}
}
\value{
If \code{plot.tree} and \code{plot.graph} are set to
\code{TRUE}, \code{extractClusters} plots two graphs. The
first shows all \eqn{\lambda}-values oredered from low to
high and indicates which values are above
\eqn{\lambda_{lim}}. Values corresponding to SOCs are
indicated in red and values corresponding to COCs are
indiced in blue. The second gives the tree illustrating
cluster merger with decreasing LD threshold where branches
represent unique clusters and branch points indicate the
merging of clusters/loci and distance gives the LD
threholds at which these events occur. Branches
corresponding to SOCs are indicated in red and values
corresponding to COCs are indiced in blue (if
\code{rm.COCs} are set to \code{FALSE}). Returns a list of
vectors giving the loci from each of the extracted
clusters.
}
\description{
to come
}
\details{
to come
}
\examples{
# Simple upper diagonal LD matrix
LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.12, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.08, 0.2, NA, NA, NA, 0.48, 0.32, 0.04, 0.44, 0.2, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.2, NA, NA, NA, NA, NA, 0.72, 0.68, 0.1, NA, NA, NA, NA, NA, NA, 0.2, 0.1, NA, NA, NA, NA, NA, NA, NA, 0.1, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
# calculate raw data
ldna <- LDnaRaw(LDmat)
# extract clusters and plot graphs, for this small example min.edges is set to zero such that each tip clade corresponds to an individual locus
par(mfcol=c(1,2))
clusters <- extractClusters(ldna, min.edges=0, phi=1)
clusters <- extractClusters(ldna, min.edges=0, phi=1, rm.COCs=FALSE)
clusters <- extractClusters(ldna, min.edges=0, phi=2, rm.COCs=FALSE)
# Bigger data set
data(LDna)
ldna <- LDnaRaw(r2.baimaii_subs)
# Different values of min.edges and phi can have crucial effects which clusters are extracted
clusters <- extractClusters(ldna, min.edges=15, phi=3)
clusters <- extractClusters(ldna, min.edges=10, phi=2)
# Includes COCs
clusters <- extractClusters(ldna, min.edges=15, phi=10, rm.COCs=FALSE)
# Only print trees
clusters <- extractClusters(ldna, min.edges=15, extract=FALSE)
clusters <- extractClusters(ldna, min.edges=5, extract=FALSE)
# Set fixed value for lambda_lim
clusters <- extractClusters(ldna, min.edges=15, lambda.lim=1)
}
\author{
Petri Kemppainen \email{petrikemppainen2@gmail.com},
Christopher Knight \email{Chris.Knight@manchester.ac.uk}
}
\seealso{
\code{\link{LDnaRaw}}
}
\keyword{extractClusters}

