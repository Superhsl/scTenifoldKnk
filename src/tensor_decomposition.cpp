#include <RcppEigen.h>
#include <vector>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;

// Compute relative Frobenius norm error between tensor and CP reconstruction
// Error = ||X - X_hat||_F / ||X||_F
double computeReconstructionError(
    const std::vector<MatrixXd>& networks,
    const MatrixXd& A,  // nGenes x K (first mode factor)
    const MatrixXd& B,  // nGenes x K (second mode factor)
    const MatrixXd& C   // nNet x K (third mode factor)
) {
    int nNet = static_cast<int>(networks.size());
    int K = static_cast<int>(C.cols());
    double totalError = 0.0;
    double totalNorm = 0.0;
    
    // Precompute rank-one matrices: A_r * B_r^T for each component
    // This avoids redundant computation in the loop
    std::vector<MatrixXd> rankOneMatrices(K);
    for (int r = 0; r < K; ++r) {
        rankOneMatrices[r].noalias() = A.col(r) * B.col(r).transpose();
    }
    
    // Compute error for each network slice
    for (int k = 0; k < nNet; ++k) {
        // Reconstruct network k: sum_r C[k,r] * A_r * B_r^T
        MatrixXd reconstructed = MatrixXd::Zero(A.rows(), B.rows());
        for (int r = 0; r < K; ++r) {
            reconstructed.noalias() += C(k, r) * rankOneMatrices[r];
        }
        
        // Accumulate squared Frobenius norms
        totalError += (networks[k] - reconstructed).squaredNorm();
        totalNorm += networks[k].squaredNorm();
    }
    
    // Return relative error
    if (totalNorm < std::numeric_limits<double>::epsilon()) {
        return 0.0;  // All-zero tensor
    }
    return std::sqrt(totalError / totalNorm);
}

// Update one factor matrix in CP-ALS (Alternating Least Squares)
// Reference: Kolda & Bader (2009). Tensor Decompositions and Applications.
//
// For mode-n update: A_n ← X_(n) * (⊙_{m≠n} A_m) * (⊙_{m≠n} A_m^T A_m)^{-1}
// where ⊙ is the Khatri-Rao product and * is the Hadamard (element-wise) product
MatrixXd updateFactor(
    const std::vector<MatrixXd>& networks,
    const MatrixXd& A,
    const MatrixXd& B,
    const MatrixXd& C,
    int mode  // 0 = update A, 1 = update B, 2 = update C
) {
    int K = static_cast<int>(A.cols());
    int nNet = static_cast<int>(networks.size());
    
    if (mode == 0) {
        // Update A (first mode factor)
        // A_new = X_(1) * (C ⊙ B) * ((C^T C) ⊙ (B^T B))^{-1}
        int nGenes = static_cast<int>(A.rows());
        MatrixXd newA = MatrixXd::Zero(nGenes, K);
        
        // Compute Gram matrices
        MatrixXd CtC = C.transpose() * C;  // K x K
        MatrixXd BtB = B.transpose() * B;  // K x K
        
        // Hadamard product of Gram matrices
        MatrixXd gram = CtC.cwiseProduct(BtB);
        gram.diagonal().array() += 1e-10;  // Regularization
        
        // Compute inverse (K is small, typically 3-5)
        MatrixXd gram_inv = gram.inverse();
        
        // Compute each row of A
        for (int i = 0; i < nGenes; ++i) {
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(K);
            
            // Sum over all network slices
            for (int k = 0; k < nNet; ++k) {
                // For each component r: rhs[r] += C[k,r] * X_k[i,:] · B[:,r]
                for (int r = 0; r < K; ++r) {
                    rhs(r) += C(k, r) * networks[k].row(i).dot(B.col(r));
                }
            }
            
            // Solve: newA[i,:] = gram^{-1} * rhs
            newA.row(i) = (gram_inv * rhs).transpose();
        }
        
        return newA;
        
    } else if (mode == 1) {
        // Update B (second mode factor)
        // B_new = X_(2) * (C ⊙ A) * ((C^T C) ⊙ (A^T A))^{-1}
        int nGenes = static_cast<int>(B.rows());
        MatrixXd newB = MatrixXd::Zero(nGenes, K);
        
        // Compute Gram matrices
        MatrixXd CtC = C.transpose() * C;
        MatrixXd AtA = A.transpose() * A;
        
        MatrixXd gram = CtC.cwiseProduct(AtA);
        gram.diagonal().array() += 1e-10;
        
        MatrixXd gram_inv = gram.inverse();
        
        for (int j = 0; j < nGenes; ++j) {
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(K);
            
            for (int k = 0; k < nNet; ++k) {
                for (int r = 0; r < K; ++r) {
                    rhs(r) += C(k, r) * networks[k].col(j).dot(A.col(r));
                }
            }
            
            newB.row(j) = (gram_inv * rhs).transpose();
        }
        
        return newB;
        
    } else {
        // Update C (third mode factor - network weights)
        // C_new[k,r] = <X_k, A[:,r] ⊗ B[:,r]> / (||A[:,r]||^2 * ||B[:,r]||^2)
        MatrixXd newC = MatrixXd::Zero(nNet, K);
        
        // Precompute denominator: ||A[:,r]||^2 * ||B[:,r]||^2 for each r
        Eigen::VectorXd denom(K);
        for (int r = 0; r < K; ++r) {
            denom(r) = A.col(r).squaredNorm() * B.col(r).squaredNorm();
        }
        denom.array() += 1e-10;  // Regularization
        
        for (int k = 0; k < nNet; ++k) {
            for (int r = 0; r < K; ++r) {
                // Compute <X_k, A[:,r] * B[:,r]^T> = B[:,r]^T * X_k^T * A[:,r]
                double inner = B.col(r).dot(networks[k].transpose() * A.col(r));
                newC(k, r) = inner / denom(r);
            }
        }
        
        return newC;
    }
}

//' CP Tensor Decomposition with C++ acceleration
//' 
//' @param networkList List of network matrices (can be sparse or dense)
//' @param K Rank of CP decomposition (number of components)
//' @param maxIter Maximum number of ALS iterations
//' @param maxError Convergence threshold (relative Frobenius norm)
//' @param nDecimal Decimal places for rounding output (0 = no rounding)
//' @return List containing:
//'   - X: reconstructed average network
//'   - A, B, C: factor matrices
//'   - iterations: number of iterations performed
//'   - error: final reconstruction error
//' @details Uses Alternating Least Squares (ALS) to compute CP decomposition
//' of the 3-way tensor formed by stacking network matrices. The reconstructed
//' network is a weighted average of the network slices.
//' @export
// [[Rcpp::export]]
List tensorDecompositionCpp(
    const List& networkList,
    int K = 3,
    int maxIter = 1000,
    double maxError = 1e-5,
    int nDecimal = 3
) {
    int nNet = static_cast<int>(networkList.size());
    
    // Input validation
    if (nNet < 1) {
        stop("networkList must contain at least 1 network");
    }
    if (K < 1) {
        stop("K (rank) must be at least 1");
    }
    if (maxIter < 1) {
        stop("maxIter must be at least 1");
    }
    if (maxError <= 0) {
        stop("maxError must be positive");
    }
    
    // Convert R list to C++ vector of dense matrices
    std::vector<MatrixXd> networks(nNet);
    
    for (int i = 0; i < nNet; ++i) {
        SEXP elem = networkList[i];
        
        if (Rf_isS4(elem)) {
            // Sparse matrix - convert to dense
            Eigen::SparseMatrix<double> sparse = as<Eigen::SparseMatrix<double>>(elem);
            networks[i] = MatrixXd(sparse);
        } else {
            networks[i] = as<MatrixXd>(elem);
        }
    }
    
    int nGenes = static_cast<int>(networks[0].rows());
    
    // Adjust K if necessary
    if (K >= nGenes) {
        K = std::max(1, nGenes - 1);
    }
    
    // Initialize factor matrices with random values using R's RNG
    // This ensures reproducibility when set.seed() is called in R
    MatrixXd A(nGenes, K);
    MatrixXd B(nGenes, K);
    MatrixXd C(nNet, K);
    
    // Use R's runif for random initialization
    for (int i = 0; i < nGenes; ++i) {
        for (int j = 0; j < K; ++j) {
            A(i, j) = R::runif(-1.0, 1.0);
            B(i, j) = R::runif(-1.0, 1.0);
        }
    }
    for (int i = 0; i < nNet; ++i) {
        for (int j = 0; j < K; ++j) {
            C(i, j) = R::runif(-1.0, 1.0);
        }
    }
    
    // Normalize columns
    for (int k = 0; k < K; ++k) {
        double normA = A.col(k).norm();
        double normB = B.col(k).norm();
        double normC = C.col(k).norm();
        
        if (normA > 1e-10) A.col(k) /= normA;
        if (normB > 1e-10) B.col(k) /= normB;
        if (normC > 1e-10) C.col(k) /= normC;
    }
    
    double prevError = 1e10;
    double currentError = 1e10;
    int iter = 0;
    
    // ALS iterations
    for (iter = 0; iter < maxIter; ++iter) {
        // Update each factor matrix in turn
        A = updateFactor(networks, A, B, C, 0);
        B = updateFactor(networks, A, B, C, 1);
        C = updateFactor(networks, A, B, C, 2);
        
        // Check convergence every 10 iterations (error computation is expensive)
        if (iter % 10 == 0) {
            currentError = computeReconstructionError(networks, A, B, C);
            
            // Convergence check
            if (currentError < maxError) {
                break;
            }
            
            // Early stopping: error reduction < 0.1%
            if (iter > 20 && (prevError - currentError) / prevError < 0.001) {
                break;
            }
            
            prevError = currentError;
        }
    }
    
    // Precompute rank-one matrices for reconstruction
    std::vector<MatrixXd> rankOneMatrices(K);
    for (int r = 0; r < K; ++r) {
        rankOneMatrices[r].noalias() = A.col(r) * B.col(r).transpose();
    }
    
    // Reconstruct weighted average network
    // CP reconstruction: X_k ≈ Σ_r C[k,r] * a_r * b_r^T
    // Average over all network slices: avg = (1/nNet) * Σ_k X_k
    MatrixXd reconstructed = MatrixXd::Zero(nGenes, nGenes);
    
    for (int k = 0; k < nNet; ++k) {
        for (int r = 0; r < K; ++r) {
            // Weight = C[k,r] / nNet (average over network slices)
            double weight = C(k, r) / static_cast<double>(nNet);
            reconstructed.noalias() += weight * rankOneMatrices[r];
        }
    }
    
    // Optional rounding (preserve precision by default)
    if (nDecimal > 0 && nDecimal <= 10) {
        double scale = std::pow(10.0, nDecimal);
        reconstructed = (reconstructed.array() * scale).round() / scale;
    }
    
    // Normalize output to [-1, 1] range
    double maxAbs = reconstructed.cwiseAbs().maxCoeff();
    if (maxAbs > 1e-10) {
        reconstructed /= maxAbs;
    }
    
    return List::create(
        Named("X") = reconstructed,
        Named("A") = A,
        Named("B") = B,
        Named("C") = C,
        Named("iterations") = iter,
        Named("error") = currentError
    );
}
