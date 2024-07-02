#' vectorSubset
#'
#' vectorSubset
#'
#' @param vec vec
#' @param mat mat
#'
#' @return matrix
#'
#' @keywords internal
vectorSubset <- function(vec, mat) {

  vmat <- c(mat)
  vvec <- vec[vmat]

  vecmat <- matrix(vvec, nrow = nrow(mat), ncol = ncol(mat))
  colnames(vecmat) <- colnames(mat)
  rownames(vecmat) <- rownames(mat)

  return(vecmat)
}
