#' @export tensorDecomposition
#' @title Performs CANDECOMP/PARAFAC (CP) Tensor Decomposition
#' @description Generate weight-averaged denoised gene regulatory networks using
#' CANDECOMP/PARAFAC (CP) Tensor Decomposition with Alternating Least Squares.
#' @param xList A list of gene regulatory networks (sparse or dense matrices)
#' @param yList Optional. A second list of gene regulatory networks for comparison
#' @param nDecimal An integer value for decimal places in output (default: 1)
#' @param K The CP rank (number of components, default: 5)
#' @param maxError Convergence threshold for relative Frobenius norm (default: 1e-5)
#' @param maxIter Maximum number of ALS iterations (default: 1000)
#' @param useCpp Logical. If TRUE (default), use faster C++ implementation
#' @return A list containing denoised network(s): $X and optionally $Y
#' @details The function stacks input networks into a 3-way tensor and performs
#' CP decomposition to extract the underlying low-rank structure. The reconstructed
#' network is a weighted average that reduces noise.
tensorDecomposition <- function(xList, yList = NULL, nDecimal = 1, K = 5,
                                maxError = 1e-5, maxIter = 1e3, useCpp = TRUE) {
  # Validate inputs
  xNets <- length(xList)
  if (xNets == 0) {
    stop("xList must contain at least one network")
  }

  # Handle yList if provided
  if (!is.null(yList)) {
    yNets <- length(yList)
    if (xNets != yNets) {
      stop("xList and yList must have the same number of networks")
    }
    nNet <- xNets
    xGenes <- unique(unlist(lapply(xList, rownames)))
    yGenes <- unique(unlist(lapply(yList, rownames)))
    sGenes <- intersect(xGenes, yGenes)
  } else {
    nNet <- xNets
    xGenes <- unique(unlist(lapply(xList, rownames)))
    sGenes <- xGenes
  }

  nGenes <- length(sGenes)
  if (nGenes == 0) {
    stop("No genes found in network matrices")
  }

  # Standardize networks to common gene set
  standardizeNetworks <- function(netList, genes) {
    lapply(netList, function(net) {
      temp <- matrix(0, length(genes), length(genes))
      rownames(temp) <- colnames(temp) <- genes
      net_mat <- as.matrix(net)
      common_genes <- genes[genes %in% rownames(net_mat)]
      temp[common_genes, common_genes] <- net_mat[common_genes, common_genes]
      Matrix::Matrix(temp, sparse = TRUE)
    })
  }

  xList_std <- standardizeNetworks(xList, sGenes)

  # Use C++ implementation if available and requested
  if (useCpp) {
    result <- tryCatch(
      {
        tensorDecompositionCpp(
          networkList = xList_std,
          K = as.integer(K),
          maxIter = as.integer(maxIter),
          maxError = maxError,
          nDecimal = as.integer(nDecimal)
        )
      },
      error = function(e) {
        warning("C++ tensor decomposition failed, falling back to R implementation: ", e$message)
        NULL
      }
    )

    if (!is.null(result)) {
      tX <- result$X
      tX <- Matrix::Matrix(tX, sparse = TRUE)
      rownames(tX) <- colnames(tX) <- sGenes

      tensorOutput <- list(X = tX)

      # Process yList if provided
      if (!is.null(yList)) {
        yList_std <- standardizeNetworks(yList, sGenes)
        resultY <- tensorDecompositionCpp(
          networkList = yList_std,
          K = as.integer(K),
          maxIter = as.integer(maxIter),
          maxError = maxError,
          nDecimal = as.integer(nDecimal)
        )
        tY <- resultY$X
        tY <- Matrix::Matrix(tY, sparse = TRUE)
        rownames(tY) <- colnames(tY) <- sGenes
        tensorOutput$Y <- tY
      }

      return(tensorOutput)
    }
  }

  # Fallback: Pure R implementation
  # Build 4D tensor (nGenes x nGenes x 1 x nNet)
  tensorX <- array(0, dim = c(nGenes, nGenes, 1, nNet))
  for (i in seq_len(nNet)) {
    tensorX[, , 1, i] <- as.matrix(xList_std[[i]])
  }

  # Run CP decomposition
  set.seed(1)
  tensorX_obj <- as.tensor(tensorX)
  decomp <- cpDecomposition(
    tnsr = tensorX_obj,
    num_components = K,
    max_iter = maxIter,
    tol = maxError
  )

  # Average reconstructed slices
  tX <- decomp$est$data[, , 1, 1]
  for (i in seq_len(nNet)[-1]) {
    tX <- tX + decomp$est$data[, , 1, i]
  }
  tX <- tX / nNet

  # Normalize and round
  max_abs <- max(abs(tX))
  if (max_abs > 0) {
    tX <- tX / max_abs
  }
  tX <- round(tX, nDecimal)
  tX <- Matrix::Matrix(tX, sparse = TRUE)
  rownames(tX) <- colnames(tX) <- sGenes

  tensorOutput <- list(X = tX)

  # Process yList if provided
  if (!is.null(yList)) {
    yList_std <- standardizeNetworks(yList, sGenes)

    tensorY <- array(0, dim = c(nGenes, nGenes, 1, nNet))
    for (i in seq_len(nNet)) {
      tensorY[, , 1, i] <- as.matrix(yList_std[[i]])
    }

    set.seed(1)
    tensorY_obj <- as.tensor(tensorY)
    decompY <- cpDecomposition(
      tnsr = tensorY_obj,
      num_components = K,
      max_iter = maxIter,
      tol = maxError
    )

    tY <- decompY$est$data[, , 1, 1]
    for (i in seq_len(nNet)[-1]) {
      tY <- tY + decompY$est$data[, , 1, i]
    }
    tY <- tY / nNet
    max_abs <- max(abs(tY))
    if (max_abs > 0) {
      tY <- tY / max_abs
    }
    tY <- round(tY, nDecimal)
    tY <- Matrix::Matrix(tY, sparse = TRUE)
    rownames(tY) <- colnames(tY) <- sGenes

    tensorOutput$Y <- tY
  }

  return(tensorOutput)
}

# Helper function to create tensor object for R implementation
as.tensor <- function(arr) {
  list(
    data = arr,
    num_modes = length(dim(arr)),
    modes = dim(arr)
  )
}
