## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 8,
  fig.height = 6,
  warning = FALSE,
  message = FALSE
)

## ----run----------------------------------------------------------------------
library(scTenifoldKnk)
library(Matrix)

data_path <- system.file("single-cell/example.csv", package = "scTenifoldKnk")
scRNAseq <- as.matrix(read.csv(data_path, row.names = 1))

result <- scTenifoldKnk(
  countMatrix = scRNAseq,
  gKO = "G100",
  qc_minLSize = 0,
  nc_nNet = 5,
  nc_nCells = 100,
  verbose = FALSE
)

## ----structure----------------------------------------------------------------
# View output structure
str(result, max.level = 2)

## ----dr_table-----------------------------------------------------------------
dr <- result$diffRegulation
knitr::kable(head(dr, 10), digits = 4,
             caption = "Top 10 Differentially Regulated Genes")

## ----distance_explain, fig.width=9, fig.height=5------------------------------
par(mfrow = c(1, 2))

# Distance distribution
hist(dr$distance, breaks = 30, col = "#3498DB", border = "white",
     main = "Distance Distribution",
     xlab = "Euclidean Distance")

# Highlight knockout gene
ko_dist <- dr$distance[dr$gene == "G100"]
abline(v = ko_dist, col = "red", lwd = 2)
text(ko_dist, par("usr")[4] * 0.9, "G100", pos = 4, col = "red")

# Distance vs Rank
plot(1:nrow(dr), dr$distance,
     pch = 19, col = "#3498DB", cex = 0.6,
     xlab = "Rank", ylab = "Distance",
     main = "Distance by Rank")
points(which(dr$gene == "G100"), ko_dist, pch = 17, col = "red", cex = 2)

## ----fc_explain, fig.width=7, fig.height=5------------------------------------
# FC distribution
par(mar = c(5, 5, 4, 2))
plot(dr$distance, dr$FC,
     pch = 19, col = "#3498DB", cex = 0.8,
     xlab = "Distance", ylab = "Fold Change",
     main = "Relationship: Distance vs Fold Change")

# The knockout gene
points(dr$distance[dr$gene == "G100"], dr$FC[dr$gene == "G100"],
       pch = 17, col = "red", cex = 2)

# Add trend
fit <- lm(FC ~ poly(distance, 2), data = dr)
x_seq <- seq(min(dr$distance), max(dr$distance), length.out = 100)
lines(x_seq, predict(fit, data.frame(distance = x_seq)), 
      col = "#E74C3C", lwd = 2)

legend("topleft", c("All genes", "G100 (KO)", "Trend"),
       pch = c(19, 17, NA), lty = c(NA, NA, 1),
       col = c("#3498DB", "red", "#E74C3C"), bty = "n")

## ----pvalue_explain, fig.width=9, fig.height=5--------------------------------
par(mfrow = c(1, 2))

# P-value vs FC
plot(dr$FC, -log10(dr$p.value + 1e-300),
     pch = 19, col = "#3498DB", cex = 0.8,
     xlab = "Fold Change", ylab = "-log10(p-value)",
     main = "FC vs Statistical Significance")
abline(h = -log10(0.05), lty = 2, col = "gray40")

# Q-Q plot for chi-square
theoretical <- qchisq(ppoints(nrow(dr)), df = 1)
observed <- sort(dr$FC)

plot(theoretical, observed,
     pch = 19, col = "#3498DB", cex = 0.6,
     xlab = expression(Theoretical~chi[1]^2),
     ylab = "Observed FC",
     main = "Chi-square Q-Q Plot")
abline(0, 1, col = "red", lwd = 2)

## ----network_compare, fig.width=10, fig.height=5------------------------------
WT <- as.matrix(result$tensorNetworks$WT)
KO <- as.matrix(result$tensorNetworks$KO)

# Focus on top affected genes
top_genes <- head(dr$gene, 10)

par(mfrow = c(1, 2))

# WT network subset
WT_sub <- WT[top_genes, top_genes]
image(1:10, 1:10, WT_sub,
      col = colorRampPalette(c("blue", "white", "red"))(100),
      axes = FALSE, main = "WT Network (Top 10 Genes)",
      xlab = "", ylab = "")
axis(1, 1:10, top_genes, las = 2, cex.axis = 0.7)
axis(2, 1:10, top_genes, las = 2, cex.axis = 0.7)

# KO network subset
KO_sub <- KO[top_genes, top_genes]
image(1:10, 1:10, KO_sub,
      col = colorRampPalette(c("blue", "white", "red"))(100),
      axes = FALSE, main = "KO Network (Top 10 Genes)",
      xlab = "", ylab = "")
axis(1, 1:10, top_genes, las = 2, cex.axis = 0.7)
axis(2, 1:10, top_genes, las = 2, cex.axis = 0.7)

## ----ko_effect, fig.width=8, fig.height=6-------------------------------------
# Compute network changes
diff_net <- WT - KO

# For knockout gene
ko_gene <- "G100"
if (ko_gene %in% rownames(WT)) {
  # Outgoing edges (regulatory targets)
  out_edges_WT <- WT[ko_gene, ]
  out_edges_KO <- KO[ko_gene, ]
  
  # Visualize
  par(mar = c(5, 5, 4, 2))
  
  plot(out_edges_WT, type = "h", col = "#3498DB", lwd = 2,
       ylim = range(c(out_edges_WT, out_edges_KO)),
       xlab = "Target Gene Index", ylab = "Edge Weight",
       main = paste("Outgoing Edges from", ko_gene))
  points(1:length(out_edges_KO), out_edges_KO, 
         type = "h", col = "#E74C3C", lwd = 2)
  
  legend("topright", c("WT", "KO"), col = c("#3498DB", "#E74C3C"),
         lwd = 2, bty = "n")
  
  cat("Sum of WT outgoing edges:", round(sum(abs(out_edges_WT)), 4), "\n")
  cat("Sum of KO outgoing edges:", round(sum(abs(out_edges_KO)), 4), "\n")
}

## ----categories, fig.width=8, fig.height=6------------------------------------
# Categorize genes
dr$category <- cut(dr$FC, 
                   breaks = c(0, 1, 2, 5, Inf),
                   labels = c("Minimal", "Moderate", "Strong", "Very Strong"))

# Visualize
par(mar = c(5, 10, 4, 2))
barplot(table(dr$category),
        col = c("#BDC3C7", "#F1C40F", "#E67E22", "#E74C3C"),
        horiz = TRUE, las = 1,
        main = "Gene Categories by Knockout Effect",
        xlab = "Number of Genes")

# Add text
cat("\nGene Categories:\n")
print(table(dr$category))

## ----qa_rank------------------------------------------------------------------
ko_rank <- which(dr$gene == "G100")
cat("Knockout gene (G100) rank:", ko_rank, "\n")

if (ko_rank == 1) {
  cat("✓ PASS: Knockout gene is top-ranked (expected behavior)\n")
} else {
  cat("⚠ WARNING: Knockout gene is not top-ranked\n")
  cat("  This may indicate:\n")
  cat("  - Weak regulatory role of the knockout gene\n")
  cat("  - Need for parameter adjustment\n")
}

## ----qa_sparsity--------------------------------------------------------------
WT_sparsity <- 1 - sum(result$tensorNetworks$WT != 0) / prod(dim(result$tensorNetworks$WT))
cat("WT network sparsity:", round(WT_sparsity, 4), "\n")

if (WT_sparsity > 0.99) {
  cat("⚠ Network may be too sparse. Consider lowering q parameter.\n")
} else if (WT_sparsity < 0.8) {
  cat("⚠ Network may be too dense. Consider increasing q parameter.\n")
} else {
  cat("✓ Network sparsity is in acceptable range.\n")
}

## ----qa_sig-------------------------------------------------------------------
n_sig_005 <- sum(dr$p.adj < 0.05)
n_sig_001 <- sum(dr$p.adj < 0.01)
pct_sig <- round(n_sig_005 / nrow(dr) * 100, 1)

cat("Significant genes (FDR < 0.05):", n_sig_005, "(", pct_sig, "%)\n")
cat("Significant genes (FDR < 0.01):", n_sig_001, "\n")

if (pct_sig > 50) {
  cat("⚠ High proportion of significant genes. May indicate strong knockout effect or need for stricter thresholds.\n")
}

## ----export, eval=FALSE-------------------------------------------------------
# # Export differential regulation table
# write.csv(dr, "diffRegulation_results.csv", row.names = FALSE)
# 
# # Export significant genes only
# sig_genes <- dr[dr$p.adj < 0.05, ]
# write.csv(sig_genes, "significant_genes.csv", row.names = FALSE)
# 
# # Export networks as edge lists
# WT_edges <- which(result$tensorNetworks$WT != 0, arr.ind = TRUE)
# edge_list <- data.frame(
#   from = rownames(result$tensorNetworks$WT)[WT_edges[, 1]],
#   to = colnames(result$tensorNetworks$WT)[WT_edges[, 2]],
#   weight = result$tensorNetworks$WT[WT_edges]
# )
# write.csv(edge_list, "WT_network_edges.csv", row.names = FALSE)

## ----session------------------------------------------------------------------
sessionInfo()

