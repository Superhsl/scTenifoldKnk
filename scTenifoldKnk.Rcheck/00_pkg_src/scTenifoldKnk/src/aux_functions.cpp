#include <RcppEigen.h>
#include "utils/matrix_ops.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;

//' Apply directional penalty to network
//' 
//' @param X Network adjacency matrix
//' @param lambda Directionality parameter (0 to 1)
//' @return Directional network matrix
//' @details For each pair (i,j), if |X[i,j]| < |X[j,i]|, then X[i,j] is set to 0.
//' The final result is: (1-lambda)*X + lambda*S where S is the directional matrix.
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> strictDirectionCpp(
    SEXP XSEXP,
    double lambda = 1.0
) {
    // Convert SEXP to Eigen matrix
    MatrixXd X = as<MatrixXd>(XSEXP);
    
    if (lambda == 0.0) {
        // Fast path: no directionality penalty, just convert to sparse
        return matrix_ops::denseToSparse(X, 1e-10);
    }
    
    // Step 1: Compute absolute values
    MatrixXd absX = X.cwiseAbs();
    MatrixXd absX_t = absX.transpose();
    
    // Step 2: Create directional mask
    // Keep X[i,j] if |X[i,j]| > |X[j,i]| (strictly greater)
    // For ties (|X[i,j]| == |X[j,i]|), keep upper triangle (i < j)
    // This matches the R implementation behavior
    int n = static_cast<int>(X.rows());
    MatrixXd mask = (absX.array() > absX_t.array()).cast<double>();
    
    // Handle ties: keep upper triangle (i < j means row index < col index)
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < j; ++i) {  // Upper triangle only
            if (std::abs(absX(i, j) - absX(j, i)) < 1e-10) {
                // Tie: keep upper triangle entry
                mask(i, j) = 1.0;
                mask(j, i) = 0.0;
            }
        }
    }
    
    // Apply mask to get directional matrix
    MatrixXd S = X.cwiseProduct(mask);
    
    // Step 3: Combine with lambda weighting
    MatrixXd result = (1.0 - lambda) * X + lambda * S;
    
    // Step 4: Convert to sparse matrix, filtering near-zero values
    return matrix_ops::denseToSparse(result, 1e-10);
}

//' Perform virtual knockout on network
//' 
//' @param network Gene regulatory network matrix
//' @param geneIdx Index of gene to knockout (0-based in C++, but R passes 1-based)
//' @return Network with gene knocked out (outgoing edges set to zero)
//' @details Sets all outgoing edges from the knocked out gene to zero,
//' meaning this gene can no longer regulate other genes.
//' @export
// [[Rcpp::export]]
Eigen::SparseMatrix<double> knockoutGeneCpp(
    SEXP networkSEXP,
    int geneIdx
) {
    // Convert SEXP to Eigen matrix
    MatrixXd network = as<MatrixXd>(networkSEXP);
    int n = static_cast<int>(network.rows());
    
    // Note: R uses 1-based indexing, C++ uses 0-based
    // The caller should handle this conversion
    if (geneIdx < 0 || geneIdx >= n) {
        stop("Invalid gene index: %d (valid range: 0 to %d)", geneIdx, n - 1);
    }
    
    // Set outgoing edges to zero (this gene cannot regulate others)
    network.row(geneIdx).setZero();
    
    // Convert to sparse
    return matrix_ops::denseToSparse(network, 1e-10);
}

//' Convert sparse matrix to dense and transpose
//' 
//' @param X Sparse matrix
//' @return Dense transposed matrix
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd sparseToDenseTranspose(
    SEXP XSEXP
) {
    // Try to convert as sparse matrix first, fallback to dense
    Eigen::MatrixXd dense;
    if (Rf_isS4(XSEXP)) {
        // Sparse matrix
        Eigen::SparseMatrix<double> sparse = as<Eigen::SparseMatrix<double>>(XSEXP);
        dense = MatrixXd(sparse);
    } else {
        // Dense matrix
        dense = as<MatrixXd>(XSEXP);
    }
    return dense.transpose();
}
