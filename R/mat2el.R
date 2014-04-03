#' Transforms matrix to an edge list.
#'
#' Transforms a lower diagonal matrix to an edge list.
#'
#' @param LDmat upper diagonal matrix of pairwise values to transform
#' @keywords mat2el
#' @export
#' @author Petri Kemppainen \email{petrikemppainen2@@gmail.com}
#' @examples
#' # Simple upper diagonal LD matrix
#' LDmat <- structure(c(NA, 0.84, 0.64, 0.24, 0.2, 0.16, 0.12, 0.44, NA, NA, 0.8, 0.28, 0.4, 0.36, 0.08, 0.4, NA, NA, NA, 0.48, 0.32, 0.04, 0.44, 0.36, NA, NA, NA, NA, 0.76, 0.56, 0.6, 0.32, NA, NA, NA, NA, NA, 0.72, 0.68, 0.28, NA, NA, NA, NA, NA, NA, 0.2, 0.24, NA, NA, NA, NA, NA, NA, NA, 0.2, NA, NA, NA, NA, NA, NA, NA, NA), .Dim = c(8L, 8L), .Dimnames = list(c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8"), c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8")))
#' # Transform
#' el <- mat2el(LDmat)
mat2el <- function(LDmat){
  from <- rownames(LDmat)[row(LDmat)[lower.tri(LDmat)]]
  to <- colnames(LDmat)[col(LDmat)[lower.tri(LDmat)]]
  weight <- LDmat[lower.tri(LDmat)]
  el <- cbind(from,to,weight)
}