pkgname <- "scTenifoldKnk"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('scTenifoldKnk')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("makeNetworksFast")
### * makeNetworksFast

flush(stderr()); flush(stdout())

### Name: makeNetworksFast
### Title: Construction of multiple gene regulatory networks
### Aliases: makeNetworksFast

### ** Examples

library(scTenifoldKnk)

# Simulating dataset
nCells <- 2000
nGenes <- 100
set.seed(1)
X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
X <- matrix(X, ncol = nCells)
rownames(X) <- paste0("gene", 1:nGenes)

# Quality control
X <- X[rowSums(X) > 0, ]

# Generate networks (sequential for small data)
networks <- makeNetworksFast(
  X = X, nNet = 10, nCells = 500, nComp = 3,
  scaleScores = TRUE, symmetric = FALSE, q = 0.95
)



cleanEx()
nameEx("scTenifoldKnk")
### * scTenifoldKnk

flush(stderr()); flush(stdout())

### Name: scTenifoldKnk
### Title: Virtual Knockout Experiment
### Aliases: scTenifoldKnk

### ** Examples

# Loading single-cell data
scRNAseq <- system.file("single-cell/example.csv", package = "scTenifoldKnk")
scRNAseq <- read.csv(scRNAseq, row.names = 1)

# Running scTenifoldKnk
result <- scTenifoldKnk(countMatrix = scRNAseq, gKO = "G100", qc_minLSize = 0)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
