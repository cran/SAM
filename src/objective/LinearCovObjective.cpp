#include <cassert>
#include "LinearCovObjective.h"

namespace SAM {

LinearCovObjective::LinearCovObjective(const double *xmat, const double *y,
                                       int n, int d, int p,
                                       bool include_intercept)
    : ObjFunction(xmat, y, n, d, p) {
  int dp = d * p;
  C.resize(dp, dp);
  Xy.resize(dp);
  Xm.resize(dp);

  // Build Gram matrix C_{jk} = X_j^T X_k / n (symmetric)
  for (int j = 0; j < d; j++) {
    for (int k = j; k < d; k++) {
      MatrixXd block = X[j].transpose() * X[k] / n;
      C.block(j * p, k * p, p, p) = block;
      if (j != k) C.block(k * p, j * p, p, p) = block.transpose();
    }
  }

  Ymean = Y.sum() / n;
  yy = Y.dot(Y) / n;

  VectorXd ones = VectorXd::Ones(n);
  for (int j = 0; j < d; j++) {
    Xy.segment(j * p, p) = X[j].transpose() * Y / n;
    Xm.segment(j * p, p) = X[j].transpose() * ones / n;
  }

  // Free X matrices — covariance update only needs C, Xy, Xm
  for (int j = 0; j < d; j++) X[j].resize(0, 0);

  if (include_intercept) model_param.intercept = Ymean;

  update_auxiliary();
  deviance = fabs(eval());
}

VectorXd LinearCovObjective::coordinate_descent(RegFunction *regfunc,
                                                int idx) {
  VectorXd beta_old = model_param.beta[idx];

  // tmp = grad[idx] + C_{idx,idx} * beta[idx]
  VectorXd tmp =
      gr[idx] + C.block(idx * p, idx * p, p, p) * model_param.beta[idx];

  model_param.beta[idx] = regfunc->threshold(tmp) * n;

  VectorXd delta = model_param.beta[idx] - beta_old;
  if (delta.squaredNorm() > 0) {
    // Covariance update: grad[k] -= C_{k,idx} * delta for all k
    for (int k = 0; k < d; k++) {
      gr[k].noalias() -= C.block(k * p, idx * p, p, p) * delta;
    }
  }
  return model_param.beta[idx];
}

void LinearCovObjective::intercept_update() {
  // mu_r = mean(r) = Ymean - sum_j Xm_j^T beta_j - intercept
  double mu_r = Ymean - model_param.intercept;
  for (int j = 0; j < d; j++)
    mu_r -= Xm.segment(j * p, p).dot(model_param.beta[j]);

  model_param.intercept += mu_r;

  // When intercept changes by mu_r, residual r -= mu_r * 1,
  // so grad[j] = X_j^T r / n changes by -Xm_j * mu_r
  for (int j = 0; j < d; j++)
    gr[j] -= Xm.segment(j * p, p) * mu_r;
}

void LinearCovObjective::update_auxiliary() {
  // Full recomputation: grad_j = Xy_j - sum_k C_{jk} beta_k - Xm_j * intercept
  for (int j = 0; j < d; j++) {
    gr[j] = Xy.segment(j * p, p) -
            Xm.segment(j * p, p) * model_param.intercept;
    for (int k = 0; k < d; k++) {
      gr[j].noalias() -= C.block(j * p, k * p, p, p) * model_param.beta[k];
    }
  }
}

void LinearCovObjective::update_gradient(int idx) {
  gr[idx] = Xy.segment(idx * p, p) -
            Xm.segment(idx * p, p) * model_param.intercept;
  for (int k = 0; k < d; k++) {
    gr[idx].noalias() -= C.block(idx * p, k * p, p, p) * model_param.beta[k];
  }
}

double LinearCovObjective::get_local_change(const VectorXd &old, int idx) {
  VectorXd tmp = old - model_param.beta[idx];
  return tmp.transpose() * C.block(idx * p, idx * p, p, p) * tmp;
}

double LinearCovObjective::get_local_change_intercept(double old) {
  double tmp = old - model_param.intercept;
  return fabs(tmp);
}

double LinearCovObjective::eval() {
  // ||y - X beta - mu||^2 / n via Gram matrix
  // = yy - 2 sum_j beta_j^T Xy_j - 2 mu Ymean
  //   + sum_{j,k} beta_j^T C_{jk} beta_k + 2 mu sum_j Xm_j^T beta_j + mu^2
  double v = yy - 2.0 * model_param.intercept * Ymean +
             model_param.intercept * model_param.intercept;
  for (int j = 0; j < d; j++) {
    v -= 2.0 * model_param.beta[j].dot(Xy.segment(j * p, p));
    v += 2.0 * model_param.intercept *
         Xm.segment(j * p, p).dot(model_param.beta[j]);
    for (int k = 0; k < d; k++) {
      v += model_param.beta[j].dot(C.block(j * p, k * p, p, p) *
                                   model_param.beta[k]);
    }
  }
  return v;
}

double LinearCovObjective::get_r2() { return n * eval(); }

}  // namespace SAM
