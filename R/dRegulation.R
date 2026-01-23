#' @export dRegulation
#' @importFrom stats dist p.adjust pchisq
#' @importFrom MASS boxcox
#' @title Differential regulation analysis
#' @description Computes differential regulation between aligned manifolds using
#' distance-based statistics. Identifies genes whose regulatory relationships
#' are most affected by the virtual knockout.
#' @param manifoldOutput Manifold alignment output matrix (2*nGenes x d dimensions)
#' @param gKO Character vector of knocked out gene name(s)
#' @return Data frame with columns: gene, distance, Z, FC, p.value, p.adj
#' @details The fold change (FC) is computed as the squared distance of each gene
#' divided by the mean squared distance of all OTHER genes (excluding the knocked
#' out gene). P-values are computed using chi-square distribution (df=1) and
#' adjusted using Benjamini-Hochberg FDR correction.
dRegulation <- function(manifoldOutput, gKO) {
  # Extract gene lists from row names
  geneList <- rownames(manifoldOutput)
  geneList <- geneList[grepl("^X_", geneList)]
  geneList <- gsub("^X_", "", geneList)
  nGenes <- length(geneList)

  eGenes <- nrow(manifoldOutput) / 2

  eGeneList <- rownames(manifoldOutput)
  eGeneList <- eGeneList[grepl("^Y_", eGeneList)]
  eGeneList <- gsub("^Y_", "", eGeneList)

  if (nGenes != eGenes) {
    stop("Number of identified and expected genes are not the same")
  }
  if (!all(eGeneList == geneList)) {
    stop("Genes are not ordered as expected. X_ genes should be followed by Y_ genes in the same order")
  }

  # Euclidean distance computation between WT and KO embeddings
  # d_i = ||X_i - Y_i||_2
  dMetric <- sqrt(rowSums((manifoldOutput[1:nGenes, , drop = FALSE] -
    manifoldOutput[(nGenes + 1):(2 * nGenes), , drop = FALSE])^2))
  names(dMetric) <- geneList

  # Box-Cox transformation for normality
  # Ensures distances are positive (required by Box-Cox)
  dMetric_positive <- pmax(dMetric, .Machine$double.eps)

  lambdaValues <- seq(-2, 2, length.out = 1000)
  lambdaValues <- lambdaValues[abs(lambdaValues) > 1e-6] # Exclude zero

  BC <- tryCatch(
    {
      result <- MASS::boxcox(dMetric_positive ~ 1, plot = FALSE, lambda = lambdaValues)
      result$x[which.max(result$y)]
    },
    error = function(e) {
      warning("Box-Cox transformation failed, using log transformation")
      0 # lambda=0 corresponds to log transformation
    }
  )

  # Apply Box-Cox transformation: (y^λ - 1) / λ
  # Note: For z-score normalization, y^λ is equivalent to (y^λ-1)/λ 
  # because they differ by only a linear transformation (scale + shift),
  # and Z = (x - mean(x))/sd(x) is invariant to linear transformations.
  # We use the power form for numerical stability.
  if (abs(BC) < 1e-6) {
    # Log transformation when lambda ~ 0 (limit of Box-Cox as λ → 0)
    nD <- log(dMetric_positive)
  } else {
    # Power transformation: y^λ (equivalent to (y^λ-1)/λ for z-score purposes)
    nD <- dMetric_positive^BC
  }

  # Z-score normalization
  Z <- as.vector(scale(nD))

  # Create output data frame
  dOut <- data.frame(
    gene = geneList,
    distance = dMetric,
    Z = Z,
    stringsAsFactors = FALSE
  )

  # Sort by distance (descending) - most affected genes first
  dOut <- dOut[order(dOut$distance, decreasing = TRUE), ]

  # Find the indices of knockout genes in the sorted data
  gKO_indices <- which(dOut$gene %in% gKO)

  if (length(gKO_indices) == 0) {
    warning(sprintf("Knockout gene(s) '%s' not found in results", paste(gKO, collapse = ", ")))
    # Use all genes for mean calculation
    reference_distances <- dOut$distance
  } else {
    # Exclude knockout gene(s) from the reference set
    reference_distances <- dOut$distance[-gKO_indices]
  }

  # Compute fold change: FC_i = d_i^2 / mean(d_j^2) for j != knockout genes
  # This measures how much gene i's distance deviates from the background
  mean_squared_dist <- mean(reference_distances^2)

  # Protect against division by zero
  if (mean_squared_dist < .Machine$double.eps) {
    warning("Mean squared distance is near zero, setting to 1 for FC calculation")
    mean_squared_dist <- 1
  }

  FC <- (dOut$distance^2) / mean_squared_dist

  # P-values using chi-square distribution with df=1
  # Under null hypothesis, FC follows chi-square(1) distribution
  pValues <- pchisq(q = FC, df = 1, lower.tail = FALSE)

  # FDR correction using Benjamini-Hochberg
  pAdjusted <- p.adjust(pValues, method = "fdr")

  dOut$FC <- FC
  dOut$p.value <- pValues
  dOut$p.adj <- pAdjusted

  return(dOut)
}
