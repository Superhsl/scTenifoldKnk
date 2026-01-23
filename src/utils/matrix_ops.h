#ifndef MATRIX_OPS_H
#define MATRIX_OPS_H

#include <RcppEigen.h>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

namespace matrix_ops {

// Compute quantile threshold using partial sorting
// Complexity: O(n) average case using nth_element
inline double computeQuantile(const Eigen::MatrixXd& X, double q) {
    // Collect non-zero absolute values
    std::vector<double> values;
    values.reserve(X.size() / 4);  // Estimate ~25% non-zero
    
    const double* data = X.data();
    const int total_size = static_cast<int>(X.size());
    
    for (int i = 0; i < total_size; ++i) {
        double val = std::abs(data[i]);
        if (val > 1e-10) {
            values.push_back(val);
        }
    }
    
    if (values.empty()) return 0.0;
    if (values.size() == 1) return values[0];
    
    // Compute quantile index
    size_t idx = static_cast<size_t>(q * (values.size() - 1));
    idx = std::min(idx, values.size() - 1);
    
    // Use nth_element for O(n) partial sorting
    std::nth_element(values.begin(), values.begin() + idx, values.end());
    
    return values[idx];
}

// Convert dense matrix to sparse with threshold filtering
inline Eigen::SparseMatrix<double> densityToSparse(
    const Eigen::MatrixXd& dense,
    double threshold
) {
    int nRows = static_cast<int>(dense.rows());
    int nCols = static_cast<int>(dense.cols());
    
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nRows * nCols / 10);  // Conservative estimate
    
    // Column-major traversal for better cache performance with Eigen's default storage
    for (int j = 0; j < nCols; ++j) {
        for (int i = 0; i < nRows; ++i) {
            double val = dense(i, j);
            if (std::abs(val) >= threshold) {
                triplets.emplace_back(i, j, val);
            }
        }
    }
    
    Eigen::SparseMatrix<double> sparse(nRows, nCols);
    sparse.setFromTriplets(triplets.begin(), triplets.end());
    sparse.makeCompressed();
    
    return sparse;
}

// Random subsample columns using reservoir sampling
// Complexity: O(nCols) time, O(nSamples) space
template<typename MatrixType>
MatrixType subsampleColumns(
    const MatrixType& X,
    int nSamples,
    unsigned int seed = 42
) {
    int nCols = static_cast<int>(X.cols());
    
    if (nSamples >= nCols) {
        return X;  // Return copy if requesting all columns
    }
    
    // Reservoir sampling algorithm (Vitter, 1985)
    std::vector<int> indices(nSamples);
    std::mt19937 gen(seed);
    
    // Fill reservoir with first nSamples indices
    for (int i = 0; i < nSamples; ++i) {
        indices[i] = i;
    }
    
    // Replace elements with decreasing probability
    for (int i = nSamples; i < nCols; ++i) {
        std::uniform_int_distribution<> dis(0, i);
        int j = dis(gen);
        if (j < nSamples) {
            indices[j] = i;
        }
    }
    
    // Sort for better cache locality during extraction
    std::sort(indices.begin(), indices.end());
    
    // Extract selected columns
    MatrixType result(X.rows(), nSamples);
    for (int i = 0; i < nSamples; ++i) {
        result.col(i) = X.col(indices[i]);
    }
    
    return result;
}

// Z-score normalization by rows (in-place)
inline void scaleRows(Eigen::MatrixXd& X) {
    int nRows = static_cast<int>(X.rows());
    
    for (int i = 0; i < nRows; ++i) {
        double mean = X.row(i).mean();
        double variance = (X.row(i).array() - mean).square().sum() / (X.cols() - 1);
        double stdDev = std::sqrt(variance);
        
        if (stdDev > 1e-10) {
            X.row(i).array() = (X.row(i).array() - mean) / stdDev;
        } else {
            X.row(i).array() = 0.0;  // Zero variance -> set to zero
        }
    }
}

// Dense to sparse conversion with threshold filtering
inline Eigen::SparseMatrix<double> denseToSparse(
    const Eigen::MatrixXd& dense,
    double threshold = 1e-10
) {
    int nRows = static_cast<int>(dense.rows());
    int nCols = static_cast<int>(dense.cols());
    
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nRows * nCols / 4);  // Estimate ~25% non-zero
    
    // Column-major traversal matches Eigen's default storage order
    for (int j = 0; j < nCols; ++j) {
        for (int i = 0; i < nRows; ++i) {
            double val = dense(i, j);
            if (std::abs(val) >= threshold) {
                triplets.emplace_back(i, j, val);
            }
        }
    }
    
    Eigen::SparseMatrix<double> sparse(nRows, nCols);
    sparse.setFromTriplets(triplets.begin(), triplets.end());
    sparse.makeCompressed();
    
    return sparse;
}

} // namespace matrix_ops

#endif // MATRIX_OPS_H
