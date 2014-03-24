\name{plotLDnetwork}
\alias{plotLDnetwork}
\title{Plots linkage disequilibrium networks}
\usage{
plotLDnetwork(ldna, el, option, threshold, clusters, summary, exl = NULL,
  full.network = TRUE, include.parent = FALSE, after.merger = FALSE)
}
\arguments{
  \item{ldna}{Output from \code{\link{LDnaRaw}}}

  \item{el}{Edge list of all pairwise LD values for
  instance produced by \code{\link{mat2el}}}

  \item{option}{\code{option=1} prints a full LD network
  from and edge list (\code{el}) and a specified LD
  threshold (\code{threshold}). \code{option=2} prints LD
  networks based on the files \code{clusters} and
  \code{summary} that contains information of extracted
  clusters.}

  \item{threshold}{Specifies the LD threshold at which an
  LD network is plotted. Only required for option=1.}

  \item{clusters}{A list of vectors giving loci for
  clusters extracted by \code{\link{extractClusters}}.}

  \item{summary}{Summary information of the clusters
  created by \code{\link{summaryLDna}} that contains.}

  \item{exl}{A list of locus names to be excluded from the
  LD networks (default is \code{NULL}).}

  \item{full.network}{If \code{TRUE} (default), includes
  all clusters in the LD networks.}

  \item{include.parent}{If \code{full.network} is set to
  \code{FALSE}, \code{include.parent} allows one to see any
  loci/clusters a given focal cluster merges with (with
  decreaisn LD threshold). When \code{FALSE}, only the
  focal cluster is shown.}

  \item{after.merger}{Whether to show LD networks at an LD
  threshold just before (\code{FALSE}) or just after
  (\code{TRUE}) merger.}
}
\description{
Allows for visual inspection of clusters extracted by
\code{\link{extractClusters}}.
}
\examples{
### Example with option=1
data(LDna)
el <- mat2el(r2.baimaii_subs)
plotLDnetwork(el=el, option=1, threshold=0.4)
### Examples with option 2
par(mfcol=c(1,2))
ldna <- LDnaRaw(r2.baimaii_subs)
clusters <- extractClusters(ldna, min.edges=20, phi=5)
summary <- summaryLDna(ldna, clusters, r2.baimaii_subs)
#default settings with option=2
plotLDnetwork(ldna, el, option=2, clusters=clusters, summary=summary)
# Examples with three useful settings (option=2)
plotLDnetwork(ldna, el, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=FALSE, after.merger=TRUE)
plotLDnetwork(ldna, el, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=TRUE, after.merger=TRUE)
plotLDnetwork(ldna, el, option=2, clusters=clusters, summary=summary, full.network=FALSE, include.parent=TRUE, after.merger=FALSE)
plotLDnetwork(ldna, el, option=2, clusters=clusters, summary=summary, full.network=TRUE, after.merger=TRUE)
}
\author{
Petri Kemppainen \email{petrikemppainen2@gmail.com}
}
\seealso{
\code{\link{LDnaRaw}}, \code{\link{extractClusters}},
\code{\link{mat2el}} and \code{\link{summaryLDna}}
}
\keyword{plotLDnetwork}

