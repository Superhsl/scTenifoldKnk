# Comprehensive Testing Suite for scTenifoldKnk v2.1.0
# Tests scientific accuracy, numerical stability, and correctness

library(scTenifoldKnk)
library(Matrix)

cat("========================================\n")
cat("  COMPREHENSIVE TESTING SUITE\n")
cat("  scTenifoldKnk v2.1.0\n")
cat("========================================\n\n")

passed <- 0
failed <- 0

test <- function(name, expr) {
  cat(sprintf("Testing: %s ... ", name))
  result <- tryCatch({
    if (isTRUE(expr)) {
      cat("PASSED\n")
      passed <<- passed + 1
      TRUE
    } else {
      cat("FAILED\n")
      failed <<- failed + 1
      FALSE
    }
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    failed <<- failed + 1
    FALSE
  })
}

# ============================================================
# 1. SCIENTIFIC ACCURACY TESTS
# ============================================================
cat("=== 1. SCIENTIFIC ACCURACY TESTS ===\n\n")

# Test 1.1: FC calculation - most critical scientific test
cat("Test 1.1: FC calculation correctness\n")
set.seed(42)
nGenes <- 10
mockManifold <- matrix(rnorm(nGenes * 2 * 2), nrow = nGenes * 2, ncol = 2)
geneNames <- paste0("G", 1:nGenes)
rownames(mockManifold) <- c(paste0("X_", geneNames), paste0("Y_", geneNames))

# Make G5 have the largest distance (knockout gene)
mockManifold[5, ] <- c(0, 0)       # X_G5
mockManifold[15, ] <- c(10, 10)    # Y_G5 - large distance

result <- dRegulation(mockManifold, "G5")

# G5 should be the top gene (highest distance)
test("G5 is top gene (highest distance)", result$gene[1] == "G5")

# Verify FC formula: FC = d^2 / mean(d^2 excluding KO)
distances <- sqrt(rowSums((mockManifold[1:nGenes, , drop = FALSE] - 
                           mockManifold[(nGenes+1):(2*nGenes), , drop = FALSE])^2))
names(distances) <- geneNames
d_G5 <- distances["G5"]
mean_other_sq <- mean(distances[names(distances) != "G5"]^2)
expected_FC_G5 <- (d_G5^2) / mean_other_sq
actual_FC_G5 <- result$FC[result$gene == "G5"]

cat(sprintf("  Expected FC for G5: %.4f\n", expected_FC_G5))
cat(sprintf("  Actual FC for G5:   %.4f\n", actual_FC_G5))
test("FC formula mathematically correct", 
     abs(expected_FC_G5 - actual_FC_G5) / expected_FC_G5 < 0.01)

cat("\n")

# Test 1.2: Distance calculation accuracy
cat("Test 1.2: Euclidean distance calculation\n")
testManifold <- matrix(0, nrow = 4, ncol = 2)
rownames(testManifold) <- c("X_A", "X_B", "Y_A", "Y_B")
testManifold["X_A", ] <- c(0, 0)
testManifold["Y_A", ] <- c(3, 4)  # Distance = 5 (3-4-5 triangle)
testManifold["X_B", ] <- c(0, 0)
testManifold["Y_B", ] <- c(1, 0)  # Distance = 1

result2 <- dRegulation(testManifold, "A")
dist_A <- result2$distance[result2$gene == "A"]
dist_B <- result2$distance[result2$gene == "B"]

cat(sprintf("  Distance A: %.6f (expected: 5)\n", dist_A))
cat(sprintf("  Distance B: %.6f (expected: 1)\n", dist_B))
test("Distance A = 5 (3-4-5 triangle)", abs(dist_A - 5) < 0.0001)
test("Distance B = 1", abs(dist_B - 1) < 0.0001)

cat("\n")

# Test 1.3: p-value validity
cat("Test 1.3: Statistical test validity\n")
test("p-values in [0,1]", all(result$p.value >= 0 & result$p.value <= 1))
test("adjusted p-values in [0,1]", all(result$p.adj >= 0 & result$p.adj <= 1))
test("FDR adjustment: p.adj >= p.value", all(result$p.adj >= result$p.value - 1e-10))

# ============================================================
# 2. NETWORK CONSTRUCTION TESTS
# ============================================================
cat("\n=== 2. NETWORK CONSTRUCTION TESTS ===\n\n")

cat("Test 2.1: pcNetFast basic functionality\n")
set.seed(123)
X <- matrix(rpois(1000, lambda = 5), nrow = 50, ncol = 20)
rownames(X) <- paste0("gene", 1:50)
X <- X[rowSums(X) > 0, ]

net <- pcNetFast(X, nComp = 3, scaleScores = TRUE, q = 0.9, verbose = FALSE)
test("Network is square", nrow(net) == ncol(net))
test("Network diagonal is zero", sum(abs(Matrix::diag(net))) < 1e-10)
test("Network is sparse dgCMatrix", inherits(net, "dgCMatrix"))
test("Gene names preserved", all(rownames(net) == rownames(X)))
test("Weights normalized to [-1,1]", max(abs(net)) <= 1 + 1e-10)

cat("\n")

cat("Test 2.2: Symmetric network option\n")
net_sym <- pcNetFast(X, nComp = 3, symmetric = TRUE, verbose = FALSE)
net_dense <- as.matrix(net_sym)
test("Symmetric: A[i,j] = A[j,i]", max(abs(net_dense - t(net_dense))) < 1e-10)

cat("\n")

cat("Test 2.3: makeNetworksFast\n")
networks <- makeNetworksFast(X, nNet = 3, nCells = 15, nComp = 3, 
                              nCores = 1, verbose = FALSE, seed = 42)
test("Correct number of networks", length(networks) == 3)
test("All networks are sparse matrices", all(sapply(networks, inherits, "dgCMatrix")))
test("All networks have correct dimensions", 
     all(sapply(networks, nrow) == nrow(X)) && all(sapply(networks, ncol) == nrow(X)))

cat("\n")

cat("Test 2.4: Reproducibility with seed\n")
nets1 <- makeNetworksFast(X, nNet = 2, nCells = 15, nComp = 3, 
                          nCores = 1, verbose = FALSE, seed = 999)
nets2 <- makeNetworksFast(X, nNet = 2, nCells = 15, nComp = 3, 
                          nCores = 1, verbose = FALSE, seed = 999)
diff1 <- sum(abs(nets1[[1]] - nets2[[1]]))
diff2 <- sum(abs(nets1[[2]] - nets2[[2]]))
cat(sprintf("  Difference net1: %.10f\n", diff1))
cat(sprintf("  Difference net2: %.10f\n", diff2))
test("Reproducible with same seed", diff1 < 1e-10 && diff2 < 1e-10)

# ============================================================
# 3. TENSOR DECOMPOSITION TESTS
# ============================================================
cat("\n=== 3. TENSOR DECOMPOSITION TESTS ===\n\n")

cat("Test 3.1: tensorDecomposition with C++\n")
td_result <- tensorDecomposition(networks, K = 2, maxIter = 50, 
                                  nDecimal = 3, useCpp = TRUE)
test("Output contains X matrix", !is.null(td_result$X))
test("X is sparse Matrix", inherits(td_result$X, "Matrix"))
test("X has correct dimensions", 
     nrow(td_result$X) == nrow(X) && ncol(td_result$X) == nrow(X))
test("X normalized to [-1,1]", max(abs(td_result$X)) <= 1 + 1e-10)

# ============================================================
# 4. MANIFOLD ALIGNMENT TESTS
# ============================================================
cat("\n=== 4. MANIFOLD ALIGNMENT TESTS ===\n\n")

cat("Test 4.1: manifoldAlignment\n")
X_net <- td_result$X
Y_net <- td_result$X  
Y_net[1, ] <- 0  # Knockout first gene

ma_result <- manifoldAlignment(X_net, Y_net, d = 2)
test("Manifold has 2*nGenes rows", nrow(ma_result) == 2 * nrow(X_net))
test("Manifold has d columns", ncol(ma_result) == 2)
test("Row names have X_/Y_ prefix", all(grepl("^X_|^Y_", rownames(ma_result))))

# ============================================================
# 5. STRICT DIRECTION TESTS
# ============================================================
cat("\n=== 5. STRICT DIRECTION TESTS ===\n\n")

cat("Test 5.1: strictDirection logic (internal function)\n")
# Using byrow=TRUE to ensure expected matrix layout
testMat <- matrix(c(0, 0.8, 0.3,
                    0.5, 0, 0.9,
                    0.2, 0.4, 0), nrow = 3, byrow = TRUE)
rownames(testMat) <- colnames(testMat) <- c("A", "B", "C")

# Access internal function using ::: 
directed <- scTenifoldKnk:::strictDirection(testMat, lambda = 1)
directed_dense <- as.matrix(directed)

# Logic: keep |X[i,j]| > |X[j,i]|, zero the other
# A->B=0.8 vs B->A=0.5: keep A->B (0.8 > 0.5)
# A->C=0.3 vs C->A=0.2: keep A->C (0.3 > 0.2)
# B->C=0.9 vs C->B=0.4: keep B->C (0.9 > 0.4)
cat(sprintf("  A->B: %.2f (should be 0.8)\n", directed_dense[1,2]))
cat(sprintf("  B->A: %.2f (should be 0)\n", directed_dense[2,1]))
cat(sprintf("  B->C: %.2f (should be 0.9)\n", directed_dense[2,3]))
cat(sprintf("  C->B: %.2f (should be 0)\n", directed_dense[3,2]))

test("Stronger edge A->B kept", directed_dense[1,2] == 0.8)
test("Weaker edge B->A zeroed", directed_dense[2,1] == 0)
test("Stronger edge B->C kept", directed_dense[2,3] == 0.9)
test("Weaker edge C->B zeroed", directed_dense[3,2] == 0)

# ============================================================
# 6. QUALITY CONTROL TESTS
# ============================================================
cat("\n=== 6. QUALITY CONTROL TESTS ===\n\n")

cat("Test 6.1: scQC basic\n")
testCounts <- matrix(rpois(500 * 100, lambda = 50), nrow = 100, ncol = 500)
rownames(testCounts) <- paste0("gene", 1:100)
colnames(testCounts) <- paste0("cell", 1:500)
rownames(testCounts)[1:5] <- paste0("MT-", LETTERS[1:5])

qc_result <- scQC(testCounts, mtThreshold = 0.2, minLSize = 100)
test("QC returns matrix", is.matrix(qc_result) || inherits(qc_result, "Matrix"))
test("QC may reduce cells", ncol(qc_result) <= ncol(testCounts))
test("QC preserves all genes", nrow(qc_result) == nrow(testCounts))

# ============================================================
# 7. INPUT VALIDATION TESTS
# ============================================================
cat("\n=== 7. INPUT VALIDATION TESTS ===\n\n")

cat("Test 7.1: pcNetFast validation\n")
badX <- matrix(1:100, ncol = 10)
rownames(badX) <- paste0("gene", 1:10)

test("Error when nComp too large", 
     tryCatch({ pcNetFast(badX, nComp = 15, verbose = FALSE); FALSE }, 
              error = function(e) grepl("nComp", e$message)))
test("Error when nComp too small",
     tryCatch({ pcNetFast(badX, nComp = 1, verbose = FALSE); FALSE }, 
              error = function(e) grepl("nComp", e$message)))

# ============================================================
# 8. END-TO-END INTEGRATION TEST
# ============================================================
cat("\n=== 8. END-TO-END INTEGRATION TEST ===\n\n")

cat("Test 8.1: Full scTenifoldKnk pipeline\n")
scRNAseq <- read.csv(system.file("single-cell/example.csv", 
                                  package = "scTenifoldKnk"), row.names = 1)
scRNAseq <- as.matrix(scRNAseq)

full_result <- tryCatch({
  scTenifoldKnk(
    countMatrix = scRNAseq,
    gKO = "G100",
    qc_minLSize = 0,
    nc_nNet = 3,
    nc_nCells = min(100, ncol(scRNAseq)),
    td_K = 3,
    td_maxIter = 100,
    verbose = FALSE
  )
}, error = function(e) {
  cat(sprintf("  Pipeline error: %s\n", e$message))
  NULL
})

test("Full pipeline completes successfully", !is.null(full_result))

if (!is.null(full_result)) {
  test("Has WT network", !is.null(full_result$tensorNetworks$WT))
  test("Has KO network", !is.null(full_result$tensorNetworks$KO))
  test("Has manifold alignment", !is.null(full_result$manifoldAlignment))
  test("Has differential regulation results", !is.null(full_result$diffRegulation))
  
  # KO gene should be most affected
  top_gene <- full_result$diffRegulation$gene[1]
  cat(sprintf("  Top affected gene: %s\n", top_gene))
  test("KO gene (G100) is most affected", top_gene == "G100")
  
  # KO network should have G100 row zeroed
  ko_net <- as.matrix(full_result$tensorNetworks$KO)
  g100_idx <- which(rownames(ko_net) == "G100")
  ko_row_sum <- sum(abs(ko_net[g100_idx, ]))
  cat(sprintf("  G100 row sum in KO network: %.10f\n", ko_row_sum))
  test("G100 outgoing edges zeroed in KO", ko_row_sum < 1e-10)
  
  # WT network should NOT have G100 row zeroed
  wt_net <- as.matrix(full_result$tensorNetworks$WT)
  wt_row_sum <- sum(abs(wt_net[g100_idx, ]))
  cat(sprintf("  G100 row sum in WT network: %.10f\n", wt_row_sum))
  test("WT network G100 row is NOT zero", wt_row_sum > 1e-10)
}

# ============================================================
# SUMMARY
# ============================================================
cat("\n========================================\n")
cat(sprintf("  RESULTS: %d PASSED, %d FAILED\n", passed, failed))
cat("========================================\n\n")

if (failed > 0) {
  cat("SOME TESTS FAILED! Please review.\n")
  quit(status = 1)
} else {
  cat("ALL TESTS PASSED! Package is scientifically accurate.\n")
}
