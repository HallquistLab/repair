// Row-wise covariance kernels for the M-MOM structural estimator.
// The dense moment algebra remains in R; these functions handle repeated
// operations over q by q by n covariance arrays.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

void check_err_cov_dimensions(
    const IntegerVector& dimensions,
    int expected_nobs,
    int expected_npredictors) {
  if (dimensions.size() != 3 ||
      dimensions[0] != expected_npredictors ||
      dimensions[1] != expected_npredictors ||
      dimensions[2] != expected_nobs) {
    stop("Internal error: err_cov must be a q by q by n array.");
  }
}

bool matrix_is_psd(const arma::mat& matrix, double tolerance) {
  arma::vec eigenvalues;
  const bool succeeded = arma::eig_sym(eigenvalues, matrix);
  if (!succeeded || eigenvalues.n_elem == 0) return false;
  const double min_eigen = eigenvalues.min();
  const double max_eigen = eigenvalues.max();
  const double scale = std::max(1.0, std::abs(max_eigen));
  return std::isfinite(min_eigen) && min_eigen >= -tolerance * scale;
}

}

// [[Rcpp::export]]
List mmom_normalize_err_cov_cpp(
    NumericVector err_cov,
    int nobs,
    int npredictors,
    double tolerance = 1e-8) {
  const IntegerVector dimensions = err_cov.attr("dim");
  const bool common_matrix = dimensions.size() == 2;
  if (common_matrix) {
    if (dimensions[0] != npredictors || dimensions[1] != npredictors) {
      stop("Internal error: common err_cov must be a q by q matrix.");
    }
  } else {
    check_err_cov_dimensions(dimensions, nobs, npredictors);
  }

  NumericVector normalized = common_matrix
      ? NumericVector(static_cast<R_xlen_t>(npredictors) * npredictors * nobs)
      : clone(err_cov);
  arma::mat row_covariance(npredictors, npredictors);
  int bad_row = 0;
  const int matrices_to_validate = common_matrix ? 1 : nobs;
  for (int row_index = 0; row_index < matrices_to_validate; ++row_index) {
    const R_xlen_t source_offset = common_matrix
        ? 0 : static_cast<R_xlen_t>(row_index) * npredictors * npredictors;
    for (int column_index = 0; column_index < npredictors; ++column_index) {
        for (
          int predictor_index = 0;
          predictor_index < npredictors;
          ++predictor_index
        ) {
        const double value = 0.5 * (
          err_cov[
            source_offset + predictor_index + npredictors * column_index
          ] +
          err_cov[
            source_offset + column_index + npredictors * predictor_index
          ]
        );
        row_covariance(predictor_index, column_index) = value;
        if (!common_matrix) {
            normalized[
              source_offset + predictor_index + npredictors * column_index
            ] = value;
        }
      }
    }
    if (!matrix_is_psd(row_covariance, tolerance)) {
      bad_row = row_index + 1;
      break;
    }
  }

  if (common_matrix && bad_row == 0) {
    for (int row_index = 0; row_index < nobs; ++row_index) {
        const R_xlen_t destination_offset =
          static_cast<R_xlen_t>(row_index) * npredictors * npredictors;
      for (int column_index = 0; column_index < npredictors; ++column_index) {
        for (
          int predictor_index = 0;
          predictor_index < npredictors;
          ++predictor_index
        ) {
            normalized[
              destination_offset + predictor_index + npredictors * column_index
            ] =
              row_covariance(predictor_index, column_index);
        }
      }
    }
  }
  normalized.attr("dim") = IntegerVector::create(
      npredictors, npredictors, nobs
  );
  return List::create(_["err_cov"] = normalized, _["bad_row"] = bad_row);
}

// [[Rcpp::export]]
int mmom_joint_err_cov_bad_row_cpp(
    const NumericVector& err_cov,
    const NumericVector& err_var_y,
    const NumericMatrix& err_cov_xy,
    double tolerance = 1e-8) {
  const IntegerVector dimensions = err_cov.attr("dim");
  if (dimensions.size() != 3) {
    stop("Internal error: err_cov must be a q by q by n array.");
  }
  const int npredictors = dimensions[0];
  const int nobs = dimensions[2];
  check_err_cov_dimensions(dimensions, nobs, npredictors);
  if (err_var_y.size() != nobs || err_cov_xy.nrow() != nobs ||
      err_cov_xy.ncol() != npredictors) {
    stop("Internal error: joint measurement-error inputs are not aligned.");
  }

  arma::mat joint(npredictors + 1, npredictors + 1, arma::fill::zeros);
  for (int row_index = 0; row_index < nobs; ++row_index) {
    joint.zeros();
    joint(0, 0) = err_var_y[row_index];
    const R_xlen_t offset =
      static_cast<R_xlen_t>(row_index) * npredictors * npredictors;
    for (int column_index = 0; column_index < npredictors; ++column_index) {
      joint(0, column_index + 1) = err_cov_xy(row_index, column_index);
      joint(column_index + 1, 0) = err_cov_xy(row_index, column_index);
        for (
          int predictor_index = 0;
          predictor_index < npredictors;
          ++predictor_index
        ) {
        joint(predictor_index + 1, column_index + 1) = err_cov[
          offset + predictor_index + npredictors * column_index
        ];
      }
    }
    if (!matrix_is_psd(joint, tolerance)) return row_index + 1;
  }
  return 0;
}

// [[Rcpp::export]]
NumericMatrix mmom_weighted_err_cov_sum_cpp(
    const NumericVector& err_cov,
    const NumericVector& weights) {
  const IntegerVector dimensions = err_cov.attr("dim");
  if (dimensions.size() != 3) {
    stop("Internal error: err_cov must be a q by q by n array.");
  }
  const int npredictors = dimensions[0];
  const int nobs = dimensions[2];
  check_err_cov_dimensions(dimensions, nobs, npredictors);
  if (weights.size() != nobs) {
    stop("Internal error: err_cov and weights are not aligned.");
  }

  NumericMatrix result(npredictors, npredictors);
  for (int row_index = 0; row_index < nobs; ++row_index) {
    const R_xlen_t offset =
      static_cast<R_xlen_t>(row_index) * npredictors * npredictors;
    for (int column_index = 0; column_index < npredictors; ++column_index) {
        for (
          int predictor_index = 0;
          predictor_index < npredictors;
          ++predictor_index
        ) {
        result(predictor_index, column_index) += weights[row_index] *
          err_cov[offset + predictor_index + npredictors * column_index];
      }
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix mmom_err_cov_matvec_cpp(
    const NumericVector& err_cov,
  const NumericVector& predictor_coef) {
  const IntegerVector dimensions = err_cov.attr("dim");
  if (dimensions.size() != 3) {
    stop("Internal error: err_cov must be a q by q by n array.");
  }
  const int npredictors = dimensions[0];
  const int nobs = dimensions[2];
  check_err_cov_dimensions(dimensions, nobs, npredictors);
  if (predictor_coef.size() != npredictors) {
    stop("Internal error: err_cov and predictor_coef are not conformable.");
  }

  NumericMatrix result(nobs, npredictors);
  for (int row_index = 0; row_index < nobs; ++row_index) {
    const R_xlen_t offset =
      static_cast<R_xlen_t>(row_index) * npredictors * npredictors;
    for (
      int predictor_index = 0;
      predictor_index < npredictors;
      ++predictor_index
    ) {
      double value = 0.0;
      for (int column_index = 0; column_index < npredictors; ++column_index) {
        value += err_cov[
          offset + predictor_index + npredictors * column_index
        ] * predictor_coef[column_index];
      }
      result(row_index, predictor_index) = value;
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericMatrix mmom_tilde_crossprod_cpp(
    const NumericVector& err_cov,
    const NumericMatrix& err_cov_xy,
    const NumericVector& coef_x,
    const NumericVector& inverse_case_variance) {
  const int nobs = err_cov_xy.nrow();
  const int npredictors = err_cov_xy.ncol();
  const IntegerVector dimensions = err_cov.attr("dim");
  check_err_cov_dimensions(dimensions, nobs, npredictors);
  if (coef_x.size() != npredictors || inverse_case_variance.size() != nobs) {
    stop("Internal error: variance-kernel inputs are not aligned.");
  }

  NumericMatrix result(npredictors + 1, npredictors + 1);
  std::vector<double> tilde(npredictors);
  for (int row_index = 0; row_index < nobs; ++row_index) {
    const R_xlen_t offset =
      static_cast<R_xlen_t>(row_index) * npredictors * npredictors;
    for (
      int predictor_index = 0;
      predictor_index < npredictors;
      ++predictor_index
    ) {
      double value = err_cov_xy(row_index, predictor_index);
      for (int column_index = 0; column_index < npredictors; ++column_index) {
        value -= err_cov[
          offset + predictor_index + npredictors * column_index
        ] * coef_x[column_index];
      }
      tilde[predictor_index] = value * inverse_case_variance[row_index];
    }
    for (int column_index = 0; column_index < npredictors; ++column_index) {
        for (
          int predictor_index = 0;
          predictor_index < npredictors;
          ++predictor_index
        ) {
        result(predictor_index + 1, column_index + 1) +=
          tilde[predictor_index] * tilde[column_index];
      }
    }
  }
  return result;
}
