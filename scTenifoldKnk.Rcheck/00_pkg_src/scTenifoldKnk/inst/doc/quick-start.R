## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----install, eval=FALSE------------------------------------------------------
# # Install from GitHub
# devtools::install_github("Zaoqu-Liu/scTenifoldKnk")

## ----load---------------------------------------------------------------------
library(scTenifoldKnk)
library(Matrix)

# Load example data
data_path <- system.file("single-cell/example.csv", package = "scTenifoldKnk")
scRNAseq <- read.csv(data_path, row.names = 1)
scRNAseq <- as.matrix(scRNAseq)

# Check data dimensions
cat("Genes:", nrow(scRNAseq), "\n")
cat("Cells:", ncol(scRNAseq), "\n")

## ----run_analysis-------------------------------------------------------------
# Perform virtual knockout of gene G100
result <- scTenifoldKnk(
  countMatrix = scRNAseq,
  gKO = "G100",
  qc_minLSize = 0,      # Skip library size filter for demo data
  nc_nNet = 5,          # Number of networks
  nc_nCells = 100,      # Cells per network
  nc_nComp = 3,         # Principal components
  td_K = 3,             # Tensor rank
  verbose = TRUE
)

## ----view_results-------------------------------------------------------------
# Top affected genes
head(result$diffRegulation, 10)

## ----visualize, fig.width=8, fig.height=6-------------------------------------
# Volcano-style plot
dr <- result$diffRegulation
dr$log10_padj <- -log10(dr$p.adj + 1e-300)
dr$significant <- dr$p.adj < 0.05

# Create plot
plot(dr$FC, dr$log10_padj,
     pch = 19,
     col = ifelse(dr$significant, "red", "gray60"),
     xlab = "Fold Change (FC)",
     ylab = "-log10(adjusted p-value)",
     main = "Virtual Knockout Effect: G100",
     cex = 0.8)

# Highlight knockout gene
ko_idx <- which(dr$gene == "G100")
points(dr$FC[ko_idx], dr$log10_padj[ko_idx], 
       pch = 17, col = "blue", cex = 2)

# Add legend
legend("topright", 
       legend = c("Significant (FDR<0.05)", "Not significant", "Knockout gene"),
       col = c("red", "gray60", "blue"),
       pch = c(19, 19, 17),
       bty = "n")

# Add threshold line
abline(h = -log10(0.05), lty = 2, col = "gray40")

## ----dr_table-----------------------------------------------------------------
knitr::kable(head(result$diffRegulation), 
             caption = "Top Differentially Regulated Genes",
             digits = 4)

## ----session------------------------------------------------------------------
sessionInfo()

