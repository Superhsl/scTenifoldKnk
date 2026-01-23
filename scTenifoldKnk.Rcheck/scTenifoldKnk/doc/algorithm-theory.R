## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  warning = FALSE,
  message = FALSE
)

## ----workflow_diagram, echo=FALSE, fig.width=10, fig.height=3-----------------
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(NULL, xlim = c(0, 10), ylim = c(0, 2), 
     xlab = "", ylab = "", axes = FALSE)

# Draw boxes
box_colors <- c("#E8F4FD", "#D4EDDA", "#FFF3CD", "#F8D7DA", "#E2D5F1", "#D1ECF1")
steps <- c("1. QC", "2. Network", "3. Tensor", "4. Knockout", "5. Manifold", "6. Diff Reg")

for (i in 1:6) {
  rect((i-1)*1.6 + 0.1, 0.5, (i-1)*1.6 + 1.5, 1.5, 
       col = box_colors[i], border = "gray40")
  text((i-1)*1.6 + 0.8, 1, steps[i], cex = 0.8, font = 2)
}

# Draw arrows
for (i in 1:5) {
  arrows((i-1)*1.6 + 1.55, 1, (i)*1.6 + 0.05, 1, 
         length = 0.1, col = "gray40", lwd = 2)
}

## ----qc_demo------------------------------------------------------------------
library(scTenifoldKnk)

# Load data
data_path <- system.file("single-cell/example.csv", package = "scTenifoldKnk")
scRNAseq <- as.matrix(read.csv(data_path, row.names = 1))

# Compute library sizes
library_sizes <- colSums(scRNAseq)

# Visualize library size distribution
hist(library_sizes, breaks = 30, col = "#4A90D9", border = "white",
     main = "Library Size Distribution",
     xlab = "Total UMI counts per cell")
abline(v = 1000, col = "red", lwd = 2, lty = 2)
legend("topright", "Threshold = 1000", col = "red", lty = 2, lwd = 2, bty = "n")

## ----pcr_illustration, fig.width=8, fig.height=4------------------------------
# Illustrate PCA on gene expression
set.seed(42)
X <- matrix(rpois(500, 10), nrow = 50, ncol = 10)
rownames(X) <- paste0("Gene", 1:50)
X <- X[rowSums(X) > 0, ]

# Standardize
X_scaled <- scale(t(X))

# PCA
pca <- prcomp(X_scaled, center = FALSE)

# Visualize
par(mfrow = c(1, 2))

# PC scores
plot(pca$x[, 1], pca$x[, 2],
     pch = 19, col = "#4A90D9",
     xlab = "PC1", ylab = "PC2",
     main = "Cell Embeddings (PC Space)")

# Variance explained
var_explained <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
barplot(var_explained[1:10], 
        names.arg = paste0("PC", 1:10),
        col = "#4A90D9",
        main = "Variance Explained",
        ylab = "% Variance",
        las = 2)

## ----tensor_demo, fig.width=8, fig.height=4-----------------------------------
# Build multiple networks
networks <- makeNetworksFast(X, nNet = 5, nCells = 8, nComp = 3, 
                              verbose = FALSE, seed = 42)

# Tensor decomposition
td_result <- tensorDecomposition(networks, K = 3, maxIter = 100)

# Visualize reconstruction
par(mfrow = c(1, 2))

# Original (first network)
net1 <- as.matrix(networks[[1]])
image(net1[1:20, 1:20], 
      col = colorRampPalette(c("blue", "white", "red"))(100),
      main = "Original Network (slice 1)",
      xlab = "Genes", ylab = "Genes", axes = FALSE)

# Reconstructed
recon <- as.matrix(td_result$X)
image(recon[1:20, 1:20],
      col = colorRampPalette(c("blue", "white", "red"))(100),
      main = "Reconstructed (Denoised)",
      xlab = "Genes", ylab = "Genes", axes = FALSE)

## ----knockout_demo, fig.width=8, fig.height=4---------------------------------
# Get WT network
WT <- as.matrix(td_result$X)

# Perform knockout
KO <- WT
knockout_gene <- "Gene1"
if (knockout_gene %in% rownames(KO)) {
  ko_idx <- which(rownames(KO) == knockout_gene)
  KO[ko_idx, ] <- 0
}

# Visualize difference
par(mfrow = c(1, 3))

# WT
image(WT[1:15, 1:15],
      col = colorRampPalette(c("blue", "white", "red"))(100),
      main = "Wild-Type Network",
      axes = FALSE)

# KO
image(KO[1:15, 1:15],
      col = colorRampPalette(c("blue", "white", "red"))(100),
      main = paste("Knockout:", knockout_gene),
      axes = FALSE)

# Difference
diff_net <- WT - KO
image(diff_net[1:15, 1:15],
      col = colorRampPalette(c("purple", "white", "orange"))(100),
      main = "Difference (WT - KO)",
      axes = FALSE)

## ----manifold_demo, fig.width=7, fig.height=6---------------------------------
# Perform manifold alignment
MA <- manifoldAlignment(Matrix::Matrix(WT, sparse = TRUE), 
                        Matrix::Matrix(KO, sparse = TRUE), 
                        d = 2)

# Extract WT and KO embeddings
n_genes <- nrow(WT)
WT_embed <- MA[1:n_genes, ]
KO_embed <- MA[(n_genes+1):(2*n_genes), ]

# Compute distances
distances <- sqrt(rowSums((WT_embed - KO_embed)^2))

# Plot
plot(WT_embed[, 1], WT_embed[, 2],
     pch = 19, col = "#4A90D9", cex = 0.8,
     xlab = "NLMA Dimension 1", ylab = "NLMA Dimension 2",
     main = "Manifold Alignment Visualization")
points(KO_embed[, 1], KO_embed[, 2],
       pch = 17, col = "#D94A4A", cex = 0.8)

# Draw connections for top affected genes
top_idx <- order(distances, decreasing = TRUE)[1:5]
for (i in top_idx) {
  lines(c(WT_embed[i, 1], KO_embed[i, 1]),
        c(WT_embed[i, 2], KO_embed[i, 2]),
        col = "orange", lwd = 2)
}

legend("topright", 
       legend = c("WT genes", "KO genes", "Top affected"),
       col = c("#4A90D9", "#D94A4A", "orange"),
       pch = c(19, 17, NA), lty = c(NA, NA, 1), lwd = c(NA, NA, 2),
       bty = "n")

## ----stats_demo, fig.width=8, fig.height=4------------------------------------
# Compute statistics
names(distances) <- rownames(WT)

# Box-Cox transformation
distances_pos <- pmax(distances, 1e-10)
bc_result <- MASS::boxcox(distances_pos ~ 1, plotit = FALSE)
lambda_opt <- bc_result$x[which.max(bc_result$y)]

# Apply transformation
if (abs(lambda_opt) < 0.01) {
  transformed <- log(distances_pos)
} else {
  transformed <- distances_pos^lambda_opt
}

# Z-scores
Z <- scale(transformed)[, 1]

par(mfrow = c(1, 2))

# Box-Cox optimization
plot(bc_result$x, bc_result$y, type = "l", lwd = 2,
     xlab = expression(lambda), ylab = "Log-Likelihood",
     main = "Box-Cox Transformation")
abline(v = lambda_opt, col = "red", lty = 2)
legend("topright", paste("Optimal λ =", round(lambda_opt, 2)),
       col = "red", lty = 2, bty = "n")

# Q-Q plot
qqnorm(Z, main = "Q-Q Plot (After Transformation)", pch = 19, col = "#4A90D9")
qqline(Z, col = "red", lwd = 2)

## ----session------------------------------------------------------------------
sessionInfo()

