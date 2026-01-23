#' @export scQC
#' @importFrom Matrix Matrix colSums
#' @importFrom stats lm predict
#' @title Single-cell quality control
#' @description Performs quality control filtering on single-cell RNA-seq data
#' based on library size, gene count, and mitochondrial content.
#' @param X A count matrix (genes x cells) or Seurat object (v4 or v5 compatible)
#' @param mtThreshold Maximum mitochondrial read ratio threshold (0-1)
#' @param minLSize Minimum library size threshold (total UMI counts per cell)
#' @return Filtered count matrix or Seurat object (same type as input)
#' @details Cells are filtered based on:
#' \enumerate{
#'   \item Library size >= minLSize
#'   \item Number of detected genes within expected range (linear model prediction)
#'   \item Mitochondrial proportion <= mtThreshold (if MT genes present)
#'   \item Library size < 2x mean (removes potential doublets)
#' }
scQC <- function(X, mtThreshold = 0.1, minLSize = 1000) {
  # Extract count matrix from Seurat object (compatible with v4 and v5)
  if (inherits(X, "Seurat")) {
    # Check Seurat version and use appropriate accessor
    countMatrix <- tryCatch(
      {
        # Seurat v5: use GetAssayData with layer argument
        if (utils::packageVersion("Seurat") >= "5.0.0") {
          SeuratObject::GetAssayData(X, layer = "counts")
        } else {
          # Seurat v4: direct slot access
          X@assays$RNA@counts
        }
      },
      error = function(e) {
        # Fallback: try GetAssayData without layer
        tryCatch(
          SeuratObject::GetAssayData(X, slot = "counts"),
          error = function(e2) {
            stop("Could not extract count matrix from Seurat object. ",
              "Please ensure Seurat and SeuratObject are properly installed.")
          }
        )
      }
    )
  } else {
    countMatrix <- X
  }

  # Ensure matrix format
  if (!inherits(countMatrix, "Matrix") && !is.matrix(countMatrix)) {
    countMatrix <- as.matrix(countMatrix)
  }

  # Step 1: Filter by minimum library size
  librarySize <- Matrix::colSums(countMatrix)
  keep_cells <- librarySize >= minLSize

  if (sum(keep_cells) == 0) {
    stop("No cells pass the minimum library size filter (minLSize = ", minLSize, ")")
  }

  countMatrix <- countMatrix[, keep_cells, drop = FALSE]
  librarySize <- Matrix::colSums(countMatrix)

  # Step 2: Compute cell-level QC metrics
  # Number of detected genes per cell
  nGenes <- Matrix::colSums(countMatrix != 0)

  # Step 3: Build linear model for expected gene count vs library size
  genesLM <- lm(nGenes ~ librarySize)
  genesLM_pred <- as.data.frame(
    predict(genesLM, data.frame(librarySize = librarySize), interval = "prediction")
  )

  # Step 4: Identify mitochondrial genes
  mtGenes <- grep("^MT-", toupper(rownames(countMatrix)), value = FALSE)

  # Step 5: Apply filtering criteria
  if (length(mtGenes) > 0) {
    # Calculate mitochondrial content
    mtCounts <- Matrix::colSums(countMatrix[mtGenes, , drop = FALSE])
    mtProportion <- mtCounts / librarySize

    # Linear model for expected MT counts
    mtLM <- lm(mtCounts ~ librarySize)
    mtLM_pred <- as.data.frame(
      predict(mtLM, data.frame(librarySize = librarySize), interval = "prediction")
    )

    # Combined filtering criteria
    selectedCells <- (mtCounts > mtLM_pred$lwr) &
      (mtCounts < mtLM_pred$upr) &
      (nGenes > genesLM_pred$lwr) &
      (nGenes < genesLM_pred$upr) &
      (mtProportion <= mtThreshold) &
      (librarySize < 2 * mean(librarySize))
  } else {
    # No MT genes found - skip MT filtering
    selectedCells <- (nGenes > genesLM_pred$lwr) &
      (nGenes < genesLM_pred$upr) &
      (librarySize < 2 * mean(librarySize))
  }

  # Get cell names that pass QC
  selectedCellNames <- colnames(countMatrix)[selectedCells]

  if (length(selectedCellNames) == 0) {
    stop("No cells pass QC filters. Consider relaxing thresholds.")
  }

  # Return filtered data in same format as input
  if (inherits(X, "Seurat")) {
    # Subset Seurat object
    X <- subset(X, cells = selectedCellNames)
  } else {
    X <- countMatrix[, selectedCellNames, drop = FALSE]
  }

  return(X)
}
