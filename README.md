# scTenifoldKnk

<!-- badges: start -->
[![R-Universe](https://zaoqu-liu.r-universe.dev/badges/scTenifoldKnk)](https://zaoqu-liu.r-universe.dev/scTenifoldKnk)
[![R-CMD-check](https://github.com/Zaoqu-Liu/scTenifoldKnk/workflows/R-CMD-check/badge.svg)](https://github.com/Zaoqu-Liu/scTenifoldKnk/actions)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
<!-- badges: end -->

## Overview

**scTenifoldKnk** is a computational framework for performing *in-silico* gene knockout experiments using single-cell RNA sequencing (scRNA-seq) data. By leveraging gene regulatory network (GRN) perturbation analysis, scTenifoldKnk enables researchers to predict the transcriptional consequences of gene knockouts without conducting wet-lab experiments.

The method integrates several advanced computational techniques:

- **Principal Component Regression (PCR)** for robust GRN inference
- **CP Tensor Decomposition** for network denoising and averaging
- **Non-linear Manifold Alignment (NLMA)** for network comparison
- **Chi-square Statistics** for identifying differentially regulated genes

## Methodology

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         scTenifoldKnk Workflow                              │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐              │
│  │ scRNA-seq│───▶│ Quality  │───▶│ Network  │───▶│  Tensor  │              │
│  │   Data   │    │ Control  │    │  Build   │    │ Decomp.  │              │
│  └──────────┘    └──────────┘    └──────────┘    └──────────┘              │
│                                                        │                    │
│                                                        ▼                    │
│  ┌──────────┐    ┌──────────┐    ┌──────────┐    ┌──────────┐              │
│  │   Diff   │◀───│ Manifold │◀───│    KO    │◀───│    WT    │              │
│  │Regulation│    │Alignment │    │ Network  │    │ Network  │              │
│  └──────────┘    └──────────┘    └──────────┘    └──────────┘              │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

The scTenifoldKnk workflow consists of six sequential steps:

1. **Quality Control**: Filter low-quality cells based on library size, gene detection rate, and mitochondrial content
2. **Network Construction**: Build multiple GRNs from random cell subsamples using principal component regression
3. **Tensor Decomposition**: Apply CANDECOMP/PARAFAC (CP) decomposition to denoise and average networks
4. **Virtual Knockout**: Simulate gene knockout by setting target gene's outgoing regulatory edges to zero
5. **Manifold Alignment**: Embed wild-type and knockout networks into a shared low-dimensional space
6. **Differential Regulation**: Identify genes significantly affected by the virtual knockout using chi-square statistics

## Installation

### From R-Universe (Recommended)

```r
install.packages("scTenifoldKnk", 
                 repos = c("https://zaoqu-liu.r-universe.dev", 
                          "https://cloud.r-project.org"))
```

### From GitHub

```r
# install.packages("devtools")
devtools::install_github("Zaoqu-Liu/scTenifoldKnk")
```

## Quick Start

```r
library(scTenifoldKnk)

# Load example dataset
data_path <- system.file("single-cell/example.csv", package = "scTenifoldKnk")
countMatrix <- as.matrix(read.csv(data_path, row.names = 1))

# Perform virtual knockout analysis
result <- scTenifoldKnk(
  countMatrix = countMatrix,
  gKO = "TargetGene",
  qc_minLSize = 1000,
  nc_nNet = 10,
  nc_nCells = 500,
  td_K = 3
)

# Examine differentially regulated genes
head(result$diffRegulation)
```

## Output Structure

The `scTenifoldKnk()` function returns a list containing:

| Component | Description |
|-----------|-------------|
| `tensorNetworks$WT` | Denoised wild-type gene regulatory network |
| `tensorNetworks$KO` | Perturbed knockout gene regulatory network |
| `manifoldAlignment` | Joint embedding matrix from NLMA |
| `diffRegulation` | Data frame of differential regulation statistics |

### Differential Regulation Table

| Column | Description |
|--------|-------------|
| `gene` | Gene identifier |
| `distance` | Euclidean distance between WT and KO embeddings |
| `Z` | Z-score (standardized distance) |
| `FC` | Fold change relative to background |
| `p.value` | Chi-square test p-value |
| `p.adj` | Benjamini-Hochberg adjusted p-value |

## Documentation

Comprehensive tutorials are available at the [package website](https://zaoqu-liu.github.io/scTenifoldKnk/):

| Vignette | Description |
|----------|-------------|
| [Quick Start](https://zaoqu-liu.github.io/scTenifoldKnk/articles/quick-start.html) | Basic usage and example workflow |
| [Algorithm Theory](https://zaoqu-liu.github.io/scTenifoldKnk/articles/algorithm-theory.html) | Mathematical foundations and methodology |
| [Visualization](https://zaoqu-liu.github.io/scTenifoldKnk/articles/visualization.html) | Result visualization techniques |
| [Advanced Usage](https://zaoqu-liu.github.io/scTenifoldKnk/articles/advanced-usage.html) | Parameter optimization and customization |
| [Result Interpretation](https://zaoqu-liu.github.io/scTenifoldKnk/articles/interpretation.html) | Biological interpretation guidelines |

## Citation

If you use scTenifoldKnk in your research, please cite:

> Osorio D, Zhong Y, Li G, Xu Q, Yang Y, Hillhouse A, Chen J, Davidson L, Tian Y, Chapkin R, Huang J, Cai JJ. **scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation.** *Patterns*. 2022;3(3):100434. doi: [10.1016/j.patter.2022.100434](https://doi.org/10.1016/j.patter.2022.100434)

## Authors

- **Zaoqu Liu** - Maintainer (<liuzaoqu@163.com>)
- **Daniel Osorio** - Original Author
- **James J. Cai** - Principal Investigator

## License

This package is licensed under GPL (>= 2).

## Related Resources

- [scTenifoldNet](https://github.com/cailab-tamu/scTenifoldNet) - Comparative single-cell GRN analysis
- [Seurat](https://satijalab.org/seurat/) - Comprehensive scRNA-seq analysis toolkit
