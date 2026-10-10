#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(scTenifoldKnk)
  library(Matrix)
})

output_dir <- Sys.getenv(
  "SCTENIFOLDKNK_OUTPUT_DIR",
  "/workspace/output"
)
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0L) stop("Refusing to overwrite non-empty output directory: ", output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

expected_versions <- c(
  scTenifoldKnk = "2.1.0",
  Matrix = "1.7.6",
  MASS = "7.3.66",
  Rcpp = "1.0.14",
  RcppEigen = "0.3.4.0.2",
  RSpectra = "0.16.2"
)
actual_versions <- vapply(
  names(expected_versions),
  function(package) as.character(packageVersion(package)),
  character(1)
)
stopifnot(identical(unname(actual_versions), unname(expected_versions)))
stopifnot(grepl("^loongarch64", R.version$platform))

data_path <- system.file("single-cell/example.csv", package = "scTenifoldKnk")
stopifnot(nzchar(data_path), file.exists(data_path))
count_matrix <- as.matrix(read.csv(data_path, row.names = 1))
stopifnot(identical(dim(count_matrix), c(100L, 3000L)))

set.seed(1)
started <- proc.time()
result <- scTenifoldKnk(
  countMatrix = count_matrix,
  gKO = "G100",
  qc_minLSize = 0,
  nc_nNet = 5,
  nc_nCells = 100,
  nc_nComp = 3,
  td_K = 3,
  verbose = TRUE
)
elapsed <- unname((proc.time() - started)[["elapsed"]])

stopifnot(
  identical(dim(result$tensorNetworks$WT), c(100L, 100L)),
  identical(dim(result$tensorNetworks$KO), c(100L, 100L)),
  identical(dim(result$manifoldAlignment), c(200L, 2L)),
  nrow(result$diffRegulation) == 100L,
  "G100" %in% result$diffRegulation$gene,
  all(c("gene", "distance", "Z", "FC", "p.value", "p.adj") %in%
        colnames(result$diffRegulation)),
  sum(result$diffRegulation$p.adj < 0.05) == 1L
)

saveRDS(result, file.path(output_dir, "result.rds"), version = 3)
write.csv(
  result$diffRegulation,
  file.path(output_dir, "diff_regulation.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session-info.txt"))
writeLines(
  c(
    "source_commit=08118337f6fc22cca5a22c443ccb78ef920dfa3b",
    paste0("package_version=", packageVersion("scTenifoldKnk")),
    paste0("r_platform=", R.version$platform),
    paste0("input_dimensions=", paste(dim(count_matrix), collapse = "x")),
    paste0("wt_dimensions=", paste(dim(result$tensorNetworks$WT), collapse = "x")),
    paste0("ko_dimensions=", paste(dim(result$tensorNetworks$KO), collapse = "x")),
    paste0("alignment_dimensions=", paste(dim(result$manifoldAlignment), collapse = "x")),
    paste0("differential_genes=", nrow(result$diffRegulation)),
    paste0("significant_fdr_0.05=", sum(result$diffRegulation$p.adj < 0.05)),
    paste0("elapsed_seconds=", format(elapsed, digits = 8)),
    "gpu_backend=none",
    "sdaa_devices=none"
  ),
  file.path(output_dir, "summary.txt")
)


cat(sprintf("Official LoongArch quick-start passed in %.3f seconds\n", elapsed))
cat("Output:", output_dir, "\n")
