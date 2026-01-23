#include <RcppEigen.h>
#include <vector>
#include <random>
#include <cmath>
#include "utils/matrix_ops.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;

// PCA using eigenvalue decomposition
// Computes top k principal components from centered data
MatrixXd fastPCA(const MatrixXd& X, int nComp) {
    int nRows = static_cast<int>(X.rows());
    int nCols = static_cast<int>(X.cols());
    
    if (nComp > std::min(nRows, nCols)) {
        nComp = std::min(nRows, nCols);
    }
    
    // Center the data (subtract column means)
    VectorXd colMeans = X.colwise().mean();
    MatrixXd Xc = X.rowwise() - colMeans.transpose();
    
    // Compute covariance matrix: Cov = X^T * X / (n-1)
    // Use noalias() to avoid temporary allocation
    MatrixXd Cov(nCols, nCols);
    Cov.noalias() = (Xc.transpose() * Xc) / (nRows - 1);
    
    // Add small regularization for numerical stability
    const double regularization = 1e-10;
    Cov.diagonal().array() += regularization;
    
    // Eigenvalue decomposition of covariance matrix
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(Cov);
    
    MatrixXd eigenvectors;
    
    if (es.info() != Eigen::Success) {
        // Fallback: use SVD (more stable but slower)
        Eigen::JacobiSVD<MatrixXd> svd(Xc, Eigen::ComputeThinV);
        eigenvectors = svd.matrixV().leftCols(nComp);
    } else {
        // Get top nComp eigenvectors (largest eigenvalues are at the END)
        eigenvectors = es.eigenvectors().rightCols(nComp);
    }
    
    // Project data onto principal components: PC = Xc * V
    MatrixXd PC(nRows, nComp);
    PC.noalias() = Xc * eigenvectors;
    
    return PC;
}

// Batch principal component regression
// Computes gene regulatory network from expression data and PC scores
// 
// Mathematical basis:
// For each gene g, we compute regression coefficients β_g = (PC^T PC)^{-1} PC^T y_g
// The network weight W[i,j] represents the shared variance between genes i and j
// through the principal component space.
MatrixXd batchPCRegression(
    const MatrixXd& Y,     // nGenes x nCells (gene expression)
    const MatrixXd& PC     // nCells x nComp (PC scores)
) {
    const int nGenes = static_cast<int>(Y.rows());
    const int nCells = static_cast<int>(Y.cols());
    const int nComp = static_cast<int>(PC.cols());
    
    if (nCells != static_cast<int>(PC.rows())) {
        stop("Dimension mismatch: Y cols (%d) != PC rows (%d)", nCells, static_cast<int>(PC.rows()));
    }
    
    // Step 1: Compute (PC^T * PC)^{-1} ONCE for all genes
    // This is an nComp x nComp matrix (typically 3x3), very small
    MatrixXd PtP(nComp, nComp);
    PtP.noalias() = PC.transpose() * PC;
    
    // Add Tikhonov regularization for numerical stability
    const double regularization = 1e-8;
    PtP.diagonal().array() += regularization;
    
    // Use Cholesky decomposition (LLT) for positive definite matrix
    Eigen::LLT<MatrixXd> llt(PtP);
    MatrixXd PtP_inv;
    
    if (llt.info() == Eigen::Success) {
        PtP_inv = llt.solve(MatrixXd::Identity(nComp, nComp));
    } else {
        // Fallback: use complete orthogonal decomposition (handles rank-deficient)
        PtP_inv = PtP.completeOrthogonalDecomposition().pseudoInverse();
    }
    
    // Step 2: Compute regression coefficients for ALL genes simultaneously
    // β = Y * PC * (PC^T PC)^{-1}
    // Dimensions: (nGenes x nCells) * (nCells x nComp) * (nComp x nComp) = (nGenes x nComp)
    MatrixXd W(nGenes, nComp);
    W.noalias() = Y * PC * PtP_inv;
    
    // Step 3: Compute network adjacency matrix
    // Network[i,j] = W[i,:] · W[j,:] = shared variance through PC space
    // Mathematically: Network = W * W^T (guaranteed symmetric)
    MatrixXd Network(nGenes, nGenes);
    Network.noalias() = W * W.transpose();
    
    return Network;
}

// Build a single network from subsampled data
MatrixXd buildSingleNetwork(
    const MatrixXd& countMatrix,
    int nCells,
    int nComp,
    double q,
    bool scaleScores,
    unsigned int seed
) {
    // Subsample cells
    MatrixXd subsample = matrix_ops::subsampleColumns(countMatrix, nCells, seed);
    
    // Transpose for PCA (we want cells as rows)
    MatrixXd X = subsample.transpose();  // nCells x nGenes
    
    // Run PCA
    MatrixXd PC = fastPCA(X, nComp);  // nCells x nComp
    
    // Batch regression to get network
    // Need to transpose back: genes as rows
    MatrixXd network = batchPCRegression(subsample, PC);
    
    // Apply quantile threshold
    double threshold = matrix_ops::computeQuantile(network, q);
    network = (network.array().abs() >= threshold).select(network, 0.0);
    
    // Scale scores if requested
    if (scaleScores) {
        double maxAbsVal = network.cwiseAbs().maxCoeff();
        if (maxAbsVal > 1e-10) {
            network /= maxAbsVal;
        }
    }
    
    return network;
}

//' Build multiple gene regulatory networks with C++ acceleration
//' 
//' @param countMatrix Gene expression matrix (genes x cells)
//' @param nNet Number of networks to generate
//' @param nCells Number of cells to subsample for each network
//' @param nComp Number of principal components
//' @param q Quantile threshold for edge filtering
//' @param scaleScores Whether to scale network weights to [-1, 1]
//' @param symmetric Whether to make networks symmetric
//' @param nThreads Number of threads (ignored, kept for API compatibility)
//' @return List of sparse network matrices
//' @details Each network is built from a random subsample of cells using
//' principal component regression. Networks are returned as sparse matrices
//' for memory efficiency.
//' @export
// [[Rcpp::export]]
List makeNetworksCpp(
    SEXP countMatrixSEXP,
    int nNet = 10,
    int nCells = 500,
    int nComp = 3,
    double q = 0.9,
    bool scaleScores = true,
    bool symmetric = false,
    int nThreads = 0  // Ignored - no OpenMP
) {
    // Convert SEXP to Eigen matrix
    MatrixXd countMatrix;
    
    if (Rf_isS4(countMatrixSEXP)) {
        // Input is sparse matrix - convert to dense for PCA
        Eigen::SparseMatrix<double> sparseMat = as<Eigen::SparseMatrix<double>>(countMatrixSEXP);
        countMatrix = MatrixXd(sparseMat);
    } else {
        countMatrix = as<MatrixXd>(countMatrixSEXP);
    }
    
    // Validate inputs
    int nGenes = static_cast<int>(countMatrix.rows());
    int totalCells = static_cast<int>(countMatrix.cols());
    
    if (nGenes < 10) {
        stop("Too few genes (< 10) for network construction");
    }
    if (totalCells < 10) {
        stop("Too few cells (< 10) for network construction");
    }
    if (nNet < 1) {
        stop("nNet must be at least 1");
    }
    if (nCells < 1) {
        stop("nCells must be at least 1");
    }
    if (nComp < 2) {
        stop("nComp must be at least 2");
    }
    if (nComp >= std::min(nGenes, nCells)) {
        nComp = std::min(nGenes, nCells) - 1;
    }
    if (q < 0.0 || q > 1.0) {
        stop("q must be between 0 and 1");
    }
    if (nCells > totalCells) {
        nCells = totalCells;
    }
    
    // Pre-allocate storage for networks (avoids repeated R memory allocation)
    std::vector<SparseMatrix<double>> networkVec(nNet);
    
    // Build networks sequentially
    for (int i = 0; i < nNet; ++i) {
        // Each network gets different random seed for reproducibility
        unsigned int seed = 1 + i;
        
        // Build single network
        MatrixXd network = buildSingleNetwork(
            countMatrix, nCells, nComp, q, scaleScores, seed
        );
        
        // Make symmetric if requested: (A + A^T) / 2
        if (symmetric) {
            network = (network + network.transpose()) * 0.5;
        }
        
        // Convert to sparse and store
        networkVec[i] = matrix_ops::densityToSparse(network, 1e-10);
    }
    
    // Convert std::vector to R List (single-threaded, safe)
    List networks(nNet);
    for (int i = 0; i < nNet; ++i) {
        networks[i] = networkVec[i];
    }
    
    return networks;
}
