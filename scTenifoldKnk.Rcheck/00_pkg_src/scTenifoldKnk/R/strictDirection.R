#' @title Apply directional penalty to network
#' @description Applies lambda penalty to enforce edge directionality.
#' For each pair of genes (i,j), keeps only the stronger edge direction.
#' @param X Network matrix (can be dense or sparse)
#' @param lambda Penalty parameter between 0 and 1. lambda=1 means full
#' directionality (keep only stronger edge), lambda=0 means no change.
#' @return Modified network matrix as sparse dgCMatrix
#' @details For each pair (i,j), if |X[i,j]| < |X[j,i]|, then X[i,j] is set to 0.
#' This enforces that regulatory relationships are directional - only the
#' stronger direction is retained.
strictDirection <- function(X, lambda = 1) {
  # Fast path: no directionality
  if (lambda == 0) {
    return(Matrix::Matrix(X, sparse = TRUE))
  }

  S <- as.matrix(X)
  n <- nrow(S)

  # Compare |X[i,j]| with |X[j,i]| for all pairs simultaneously
  absS <- abs(S)
  absS_t <- t(absS)

  # Create mask: TRUE where |X[i,j]| >= |X[j,i]|
  # This means we keep X[i,j], zero out X[j,i]
  keep_mask <- absS >= absS_t

  # For ties (|X[i,j]| == |X[j,i]|), keep upper triangle convention
  # This ensures we don't zero out both directions for equal values
  tie_mask <- (absS == absS_t)
  upper_tri <- upper.tri(S, diag = FALSE)
  keep_mask[tie_mask & upper_tri] <- TRUE
  keep_mask[tie_mask & !upper_tri & t(tie_mask & upper_tri)] <- FALSE

  # Apply mask to create directional matrix
  S_directed <- S * keep_mask

  # Combine with original based on lambda
  O <- (1 - lambda) * S + lambda * S_directed

  # Convert to sparse matrix
  O <- Matrix::Matrix(O, sparse = TRUE)

  return(O)
}
