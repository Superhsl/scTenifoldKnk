#include <RcppEigen.h>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <numeric>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Box-Cox transformation with COMPLETE profile log-likelihood
// Reference: Box, G. E., & Cox, D. R. (1964). An analysis of transformations.
// Full log-likelihood: LL(λ) = -n/2 * log(σ²) + (λ-1) * Σlog(y_i)
VectorXd boxCoxTransform(const VectorXd& x) {
    int n = static_cast<int>(x.size());
    
    // Critical: Box-Cox transformation requires strictly positive values
    double minVal = x.minCoeff();
    VectorXd xPositive = x;
    
    if (minVal <= 0) {
        // Add offset to make all values positive
        double offset = std::abs(minVal) + 1e-6;
        xPositive = x.array() + offset;
    }
    
    // Additional check: ensure no values are too close to zero
    xPositive = xPositive.array().max(1e-10);
    
    // Precompute sum of log(y) for Jacobian term (constant across lambda)
    double sumLogY = xPositive.array().log().sum();
    
    // Find optimal lambda using COMPLETE profile log-likelihood
    // Search range: -2 to 2, excluding values near 0
    const int nLambda = 1000;
    VectorXd lambdaValues = VectorXd::LinSpaced(nLambda, -2.0, 2.0);
    VectorXd logLikelihood(nLambda);
    
    for (int i = 0; i < nLambda; ++i) {
        double lambda = lambdaValues(i);
        
        VectorXd transformed(n);
        
        if (std::abs(lambda) < 1e-6) {
            // Log transformation when lambda ≈ 0
            // Box-Cox: y^(λ) = (y^λ - 1)/λ → log(y) as λ → 0
            transformed = xPositive.array().log();
        } else {
            // Standard Box-Cox transformation: (y^λ - 1)/λ
            // For likelihood, we use y^λ (equivalent up to linear transformation)
            transformed = xPositive.array().pow(lambda);
        }
        
        // Compute variance of transformed data
        double mean = transformed.mean();
        double variance = (transformed.array() - mean).square().sum() / (n - 1);
        
        // Protect against zero variance
        if (variance < 1e-20) variance = 1e-20;
        
        // COMPLETE log-likelihood with Jacobian term
        // LL(λ) = -n/2 * log(σ²) + (λ-1) * Σlog(y_i)
        double ll = -0.5 * n * std::log(variance) + (lambda - 1.0) * sumLogY;
        logLikelihood(i) = ll;
    }
    
    // Find lambda with maximum likelihood
    int maxIdx;
    logLikelihood.maxCoeff(&maxIdx);
    double optimalLambda = lambdaValues(maxIdx);
    
    // Apply optimal transformation
    VectorXd transformed(n);
    if (std::abs(optimalLambda) < 1e-6) {
        // Log transformation
        transformed = xPositive.array().log();
    } else {
        transformed = xPositive.array().pow(optimalLambda);
    }
    
    // Validate output: check for NaN or Inf
    bool hasInvalid = false;
    for (int i = 0; i < n; ++i) {
        if (std::isnan(transformed(i)) || std::isinf(transformed(i))) {
            hasInvalid = true;
            break;
        }
    }
    
    if (hasInvalid) {
        // Fallback to log transformation
        return xPositive.array().log();
    }
    
    return transformed;
}

// Compute chi-square p-values
VectorXd chiSquarePvalues(const VectorXd& chiSquareStats, int df = 1) {
    int n = static_cast<int>(chiSquareStats.size());
    VectorXd pvalues(n);
    
    for (int i = 0; i < n; ++i) {
        // P[X > x] for chi-square distribution
        // Using R's pchisq via Rcpp (lower.tail = FALSE, log.p = FALSE)
        pvalues(i) = R::pchisq(chiSquareStats(i), static_cast<double>(df), 0, 0);
    }
    
    return pvalues;
}

// FDR correction (Benjamini-Hochberg)
VectorXd fdrCorrection(const VectorXd& pvalues) {
    int n = static_cast<int>(pvalues.size());
    
    // Create index vector
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    
    // Sort indices by p-values
    std::sort(indices.begin(), indices.end(), [&](int a, int b) {
        return pvalues(a) < pvalues(b);
    });
    
    VectorXd adjusted(n);
    double cumMin = 1.0;
    
    // BH procedure (from largest to smallest)
    for (int i = n - 1; i >= 0; --i) {
        int idx = indices[i];
        double bh = pvalues(idx) * n / (i + 1);
        cumMin = std::min(cumMin, bh);
        adjusted(idx) = std::min(cumMin, 1.0);
    }
    
    return adjusted;
}

//' Differential regulation analysis with C++ acceleration
//' 
//' @param manifoldOutput Matrix with WT and KO gene embeddings (2*nGenes x d)
//' @param geneNames Vector of gene names
//' @param gKO Name of knocked out gene
//' @return DataFrame with differential regulation statistics
//' @export
// [[Rcpp::export]]
DataFrame dRegulationCpp(
    SEXP manifoldOutputSEXP,
    const std::vector<std::string>& geneNames,
    const std::string& gKO
) {
    // Convert SEXP to Eigen matrix
    MatrixXd manifoldOutput = as<MatrixXd>(manifoldOutputSEXP);
    
    int nGenes = static_cast<int>(geneNames.size());
    
    // Validate dimensions
    if (static_cast<int>(manifoldOutput.rows()) != 2 * nGenes) {
        stop("manifoldOutput should have 2*nGenes rows");
    }
    
    // Compute Euclidean distances between WT and KO embeddings
    // WT embeddings: rows 0 to nGenes-1
    // KO embeddings: rows nGenes to 2*nGenes-1
    MatrixXd wtEmbeddings = manifoldOutput.topRows(nGenes);
    MatrixXd koEmbeddings = manifoldOutput.bottomRows(nGenes);
    
    // Distance: ||WT_i - KO_i||_2 for each gene
    VectorXd distances = (wtEmbeddings - koEmbeddings).rowwise().norm();
    
    // Box-Cox transformation for normality
    VectorXd transformed = boxCoxTransform(distances);
    
    // Z-score normalization
    double mean = transformed.mean();
    double variance = (transformed.array() - mean).square().sum() / (nGenes - 1);
    double stdDev = std::sqrt(variance);
    
    // Prevent division by zero
    if (stdDev < std::numeric_limits<double>::epsilon()) {
        stdDev = 1.0;
    }
    
    VectorXd Z = (transformed.array() - mean) / stdDev;
    
    // Sort genes by distance (descending) - most affected genes first
    std::vector<int> sortedIdx(nGenes);
    std::iota(sortedIdx.begin(), sortedIdx.end(), 0);
    std::sort(sortedIdx.begin(), sortedIdx.end(), [&](int a, int b) {
        return distances(a) > distances(b);
    });
    
    // Find ALL knockout gene indices (support multiple KO genes)
    std::vector<int> koIndices;
    for (int i = 0; i < nGenes; ++i) {
        if (geneNames[i] == gKO) {
            koIndices.push_back(i);
        }
    }
    
    // Compute FC: distance^2 / mean(reference distances^2)
    // Reference = all genes EXCEPT the knockout gene(s)
    VectorXd distSquared = distances.array().square();
    
    // Compute mean squared distance excluding KO gene(s)
    double sumSquared = distSquared.sum();
    int referenceCount = nGenes;
    
    for (int idx : koIndices) {
        sumSquared -= distSquared(idx);
        referenceCount--;
    }
    
    // Handle edge case: all genes are KO genes
    if (referenceCount <= 0) {
        referenceCount = nGenes;
        sumSquared = distSquared.sum();
    }
    
    double meanSquared = sumSquared / referenceCount;
    
    // Protect against division by zero
    if (meanSquared < std::numeric_limits<double>::epsilon()) {
        meanSquared = 1.0;
    }
    
    VectorXd FC = distSquared / meanSquared;
    
    // P-values using chi-square distribution (df=1)
    VectorXd pvalues = chiSquarePvalues(FC, 1);
    
    // FDR correction using Benjamini-Hochberg
    VectorXd padj = fdrCorrection(pvalues);
    
    // Create output DataFrame (sorted by distance, descending)
    std::vector<std::string> sortedGenes(nGenes);
    VectorXd sortedDist(nGenes);
    VectorXd sortedZ(nGenes);
    VectorXd sortedFC(nGenes);
    VectorXd sortedPval(nGenes);
    VectorXd sortedPadj(nGenes);
    
    for (int i = 0; i < nGenes; ++i) {
        int idx = sortedIdx[i];
        sortedGenes[i] = geneNames[idx];
        sortedDist(i) = distances(idx);
        sortedZ(i) = Z(idx);
        sortedFC(i) = FC(idx);
        sortedPval(i) = pvalues(idx);
        sortedPadj(i) = padj(idx);
    }
    
    DataFrame result = DataFrame::create(
        Named("gene") = sortedGenes,
        Named("distance") = sortedDist,
        Named("Z") = sortedZ,
        Named("FC") = sortedFC,
        Named("p.value") = sortedPval,
        Named("p.adj") = sortedPadj
    );
    
    return result;
}
