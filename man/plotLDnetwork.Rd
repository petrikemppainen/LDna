\name{plotLDnetwork}
\alias{plotLDnetwork}
\title{Plots linkage disequilibrium networks}
\usage{
plotLDnetwork(ldna, LDmat, option, threshold, clusters, summary, exl = NULL,
  full.network = TRUE, include.parent = FALSE, after.merger = FALSE,
  graph.object = FALSE)
}
\arguments{
  \item{ldna}{Output from \code{\link{LDnaRaw}}}

  \item{LDmat}{Lower diagonal LD matrix}

  \item{option}{\code{option=1} prints a full LD network
  from and edge list (\code{LDmat}) and a specified LD
  threshold (\code{threshold}). \code{option=2} prints LD
  networks based on the files \code{clusters} and
  \code{summary} that contains information of extracted
  clusters.}

  \item{threshold}{Specifies the LD threshold at which an
  LD network is plotted. Only required for option=1.}

  \item{clusters}{Output from
  \code{\link{extractClusters}}.}

  \item{summary}{Output from \code{\link{summaryLDna}}.}

  \item{exl}{A list of locus names to be excluded from the
  LD networks (default is \code{NULL}).}

  \item{full.network}{If \code{TRUE} (default), includes
  all loci in the LD networks, not recommended for large
  data sets.}

  \item{include.parent}{If \code{full.network=FALSE} and
  \code{include.parent=TRUE} all loci from the parent
  cluster (after merger) are included. If
  \code{include.parent=FALSE} only the focal cluster is
  shown.}

  \item{after.merger}{Whether to show LD networks at an LD
  threshold just before (\code{FALSE}) or just after
  (\code{TRUE}) merger.}

  \item{graph.object}{Whether to output \code{graph object}
  when option=1 (default=\code{FALSE}).}
}
\description{
Allows for visual inspection of clusters identified by
\code{\link{extractClusters}}.
}
\examples{
### Example with option=1
data(LDna)
plotLDnetwork(LDmat=r2.baimaii_subs, option=1, threshold=0.4)
### Examples with option 2
par(mfcol=c(1,2))
ldna <- LDnaRaw(r2.baimaii_subs)
clusters <- extractClusters(ldna, min.edges=20, phi=5)
summary <- summaryLDna(ldna, clusters, r2.baimaii_subs)
#default settings with option=2
plotLDnetwork(ldna, r2.baimaii_subs, option=2, clusters=clusters, summary=summary)
## Other useful settings
# For large data sets
plotLDnetwork(ldna, r2.baimaii_subs, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=FALSE, after.merger=FALSE)
# To visualise the merger
plotLDnetwork(ldna, r2.baimaii_subs, option=2, clusters=clusters, summary=summary, full.network=TRUE, after.merger=TRUE)
# Or
plotLDnetwork(ldna, r2.baimaii_subs, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=TRUE, after.merger=TRUE)
# To show that ususally several clusters are involved in most mergers
plotLDnetwork(ldna, r2.baimaii_subs, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=TRUE, after.merger=FALSE)
}
\author{
Petri Kemppainen \email{petrikemppainen2@gmail.com}
}
\seealso{
\code{\link{LDnaRaw}}, \code{\link{extractClusters}} and
\code{\link{summaryLDna}}
}
\keyword{plotLDnetwork}

