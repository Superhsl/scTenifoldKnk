#' @export makeNetworksFast
#' @importFrom Matrix Matrix rowSums
#' @importFrom parallel detectCores makeCluster clusterExport parLapply stopCluster clusterEvalQ
#' @title Construction of multiple gene regulatory networks
#' @description Computes multiple gene regulatory networks from random subsamples of cells using
#' principal component regression. Supports both sequential and parallel processing.
#' @param X A filtered and normalized gene expression matrix with cells as columns and genes as rows.
#' @param nNet An integer value. The number of networks to generate (default: 10).
#' @param nCells An integer value. The number of cells to subsample for each network (default: 500).
#' @param nComp An integer value. The number of principal components (default: 3).
#' Must be >= 2 and < number of genes.
#' @param scaleScores Logical. If TRUE, normalize weights so max absolute value is 1.
#' @param symmetric Logical. If TRUE, return symmetric network matrices.
#' @param q A decimal value between 0 and 1. Quantile threshold for edge filtering.
#' @param nCores An integer value. Number of cores for parallel processing. Set to 1 for sequential.
#' Default: detectCores() - 1.
#' @param verbose Logical. If TRUE, print progress information.
#' @param seed An integer value. Random seed for reproducibility. If NULL, no seed is set.
#' @return A list of nNet gene regulatory networks in dgCMatrix format.
#' @details Each network is constructed independently from a random subsample of cells,
#' making the process embarrassingly parallel. The function automatically chooses
#' sequential processing for small datasets where parallel overhead exceeds benefits.
#' @examples
#' library(scTenifoldKnk)
#'
#' # Simulating dataset
#' nCells <- 2000
#' nGenes <- 100
#' set.seed(1)
#' X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
#' X <- matrix(X, ncol = nCells)
#' rownames(X) <- paste0("gene", 1:nGenes)
#'
#' # Quality control
#' X <- X[rowSums(X) > 0, ]
#'
#' # Generate networks (sequential for small data)
#' networks <- makeNetworksFast(
#'   X = X, nNet = 10, nCells = 500, nComp = 3,
#'   scaleScores = TRUE, symmetric = FALSE, q = 0.95
#' )
makeNetworksFast <- function(X,
                             nNet = 10,
                             nCells = 500,
                             nComp = 3,
                             scaleScores = TRUE,
                             symmetric = FALSE,
                             q = 0.95,
                             nCores = NULL,
                             verbose = TRUE,
                             seed = 1) {
  # Input validation and coercion
  if (!is.matrix(X) && !inherits(X, "dgCMatrix")) {
    X <- as.matrix(X)
  }

  geneList <- rownames(X)
  if (is.null(geneList) || length(geneList) == 0) {
    geneList <- paste0("gene", seq_len(nrow(X)))
    rownames(X) <- geneList
  }
  nGenes <- nrow(X)
  nCol <- ncol(X)

  # Validate parameters
  if (nGenes == 0) {
    stop("No genes in input matrix")
  }
  if (nCells > nCol) {
    warning(sprintf(
      "nCells (%d) > available cells (%d), using all cells",
      nCells, nCol
    ))
    nCells <- nCol
  }
  if (nComp >= nGenes) {
    stop("nComp must be less than the number of genes")
  }
  if (nComp < 2) {
    stop("nComp must be >= 2")
  }

  # Determine number of cores
  if (is.null(nCores)) {
    nCores <- max(1, parallel::detectCores() - 1)
  }
  nCores <- min(nCores, nNet)  # No point having more cores than networks

  # Smart parallelization decision
  # Parallel overhead ~1-2s; only worthwhile for larger jobs
  use_parallel <- (nGenes > 200 || nNet > 20) && nCores > 1

  if (!use_parallel && nCores > 1) {
    if (verbose) {
      cat("Using sequential processing (dataset size doesn't warrant parallel overhead)\n")
    }
    nCores <- 1
  }

  if (verbose) {
    cat(sprintf(
      "Generating %d networks (%d genes, %d cells/network, %d PCs)...\n",
      nNet, nGenes, nCells, nComp
    ))
  }

  # Set seed for reproducibility BEFORE generating samples
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Pre-generate all random samples (ensures reproducibility)
  cell_samples <- lapply(seq_len(nNet), function(i) {
    sample(x = seq_len(nCol), size = nCells, replace = TRUE)
  })

  # Convert to dense matrix for processing
  X_matrix <- as.matrix(X)

  # Network builder function (defined once, used in both paths)
  buildOneNetwork <- function(net_idx) {
    selected_cells <- cell_samples[[net_idx]]
    X_sub <- X_matrix[, selected_cells, drop = FALSE]

    # Filter zero-expression genes in this subsample
    gene_sums <- rowSums(X_sub)
    X_sub <- X_sub[gene_sums > 0, , drop = FALSE]

    if (nrow(X_sub) < nComp + 1) {
      warning(sprintf(
        "Network %d: too few genes (%d) for nComp=%d, skipping",
        net_idx, nrow(X_sub), nComp
      ))
      return(NULL)
    }

    # Build network
    net <- pcNetFast(
      X_sub,
      nComp = nComp,
      scaleScores = scaleScores,
      symmetric = symmetric,
      q = q,
      verbose = FALSE
    )

    # Expand to full gene set
    O <- matrix(0, nrow = nGenes, ncol = nGenes)
    rownames(O) <- colnames(O) <- geneList
    active_genes <- rownames(net)
    O[active_genes, active_genes] <- as.matrix(net)
    O <- Matrix::Matrix(O, sparse = TRUE)

    return(O)
  }

  # Build networks
  if (nCores > 1 && nNet > 1) {
    if (verbose) cat(sprintf("Building networks in parallel (%d cores)...\n", nCores))

    # Set up cluster with proper cleanup on exit
    cl <- parallel::makeCluster(nCores)
    on.exit(parallel::stopCluster(cl), add = TRUE)  # Ensures cleanup

    # Export to workers
    parallel::clusterExport(
      cl,
      varlist = c(
        "X_matrix", "geneList", "nGenes", "nComp",
        "scaleScores", "symmetric", "q", "cell_samples"
      ),
      envir = environment()
    )

    # Load packages on workers
    parallel::clusterEvalQ(cl, {
      library(Matrix)
      library(RSpectra)
      library(scTenifoldKnk)
    })

    # Parallel build
    networks <- parallel::parLapply(cl, seq_len(nNet), buildOneNetwork)

  } else {
    if (verbose) cat("Building networks sequentially...\n")

    networks <- lapply(seq_len(nNet), function(i) {
      if (verbose && i %% 5 == 1) {
        cat(sprintf("  Network %d/%d\n", i, nNet))
      }
      buildOneNetwork(i)
    })
  }

  # Remove failed networks
  null_mask <- sapply(networks, is.null)
  if (any(null_mask)) {
    warning(sprintf("%d/%d networks failed to build", sum(null_mask), nNet))
    networks <- networks[!null_mask]
  }

  if (length(networks) == 0) {
    stop("All networks failed to build. Check input data and parameters.")
  }

  if (verbose) {
    cat(sprintf("Successfully generated %d networks\n", length(networks)))
  }

  return(networks)
}
