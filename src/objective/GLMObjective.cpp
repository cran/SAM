#include "GLMObjective.h"
#include "../utils.h"
#include <stdio.h>
#include <iostream>
#include <Rcpp.h>
#include <RcppEigen.h>
#include <complex>
#include <algorithm>

//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

namespace SAM {
  const double eps = 1e-4;
  using Eigen::ArrayXd;
  using Eigen::VectorXcd;

  GLMObjective::GLMObjective(const double *xmat, const double *y, int n, int d, int p,
                             double step_size0, bool include_intercept)
    : ObjFunction(xmat, y, n, d, p), P(n), W(n), R(n), sum_r(0), sum_w(0), step_size0(step_size0),
      wxx_gen(0), wxx_cached_gen(d, -1), wxx_cache(d), step_cache(d, 0.0) {

    if (include_intercept) {
      double avr_y = Y.sum() / n;

      model_param.intercept = log(avr_y / (1 - avr_y));
    }
  }


  VectorXd GLMObjective::coordinate_descent(RegFunction *regfunc, int idx) {
    VectorXd gr = X[idx].transpose() * -R / n;
    VectorXd tmp;

    MatrixXd wXX;
    double step_size;

    if (wxx_cached_gen[idx] == wxx_gen) {
      // Reuse cached wXX and step_size (W unchanged since last computation)
      wXX = wxx_cache[idx];
      step_size = step_cache[idx];
    } else {
      wXX.resize(p, p);
      wXX.setZero();
      for (int i = 0; i < n; i++)
        wXX += W[i] * X[idx].row(i).transpose() * X[idx].row(i);
      wXX /= n;

      VectorXcd eigenvalues = wXX.eigenvalues();
      step_size = 0;
      for (int i = 0; i < eigenvalues.size(); i++) {
        step_size = std::max(step_size, eigenvalues[i].real());
      }

      wxx_cache[idx] = wXX;
      step_cache[idx] = step_size;
      wxx_cached_gen[idx] = wxx_gen;
    }
    assert(step_size >= 0);

    tmp = regfunc->threshold_p((model_param.beta[idx] - gr / step_size), step_size);
    VectorXd delta_beta = tmp - model_param.beta[idx];
    tmp = (wXX*model_param.beta[idx] - gr);
    VectorXd old_beta = model_param.beta[idx];
    model_param.beta[idx] = regfunc->threshold(tmp)/step_size;

    delta_beta = model_param.beta[idx] - old_beta;

    if (calc_norm(delta_beta) > 1e-8) {
      Xb += X[idx] * delta_beta;

      R -= W.cwiseProduct(X[idx] * delta_beta);
    }

    return model_param.beta[idx];
  }

  void GLMObjective::intercept_update() {
    sum_r = R.sum();
    model_param.intercept += sum_r/sum_w;
    R -= sum_r/sum_w * W;
    sum_r = 0;
  }

  void GLMObjective::update_gradient(int idx) {
    gr[idx] = X[idx].transpose() * (P - Y) / n;
  }

  double GLMObjective::get_local_change(const VectorXd &old, int idx) {
    VectorXd delta_beta = old - model_param.beta[idx];

    double max_eigen;
    if (wxx_cached_gen[idx] == wxx_gen) {
      // step_cache = max_eig(wXX/n), need max_eig(wXX) = step_cache * n
      max_eigen = step_cache[idx] * n;
    } else {
      MatrixXd wXX(p, p);
      wXX.setZero();
      for (int i = 0; i < n; i++)
        wXX += W[i] * X[idx].row(i).transpose() * X[idx].row(i);
      VectorXcd eigenvalues = wXX.eigenvalues();
      max_eigen = 0;
      for (int i = 0; i < eigenvalues.size(); i++) {
        max_eigen = std::max(max_eigen, eigenvalues[i].real());
      }
    }
    return max_eigen * (delta_beta.dot(delta_beta)) / (2 * n);
  }
  double GLMObjective::get_local_change_intercept(double old) {
    double tmp = old - model_param.intercept;
    return (sum_w * tmp * tmp / (2 * n));
  }
  double GLMObjective::get_r2() {
    // GLM families use deviance for model fit; sse[] output is not meaningful.
    return 0;
  }

  LogisticObjective::LogisticObjective(const double *xmat, const double *y, int n,
                                       int d, int p, double step_size0, bool include_intercept)
    : GLMObjective(xmat, y, n, d, p, step_size0, include_intercept) {
    update_auxiliary();
    for (int i = 0; i < d; i++) update_gradient(i);

    model_param.intercept = 0.0;
    update_auxiliary();

    deviance = fabs(eval());
  }

  void LogisticObjective::update_auxiliary() {
    P = -(Xb.array() + model_param.intercept);
    P = P.array().exp();
    P = (P.array() + 1.0).inverse();
    R = Y - P;

    W = P.array() * -(P.array() - 1);
    sum_w = W.sum();
    wxx_gen++;
  }

  double LogisticObjective::eval() {
    double v = 0.0;

    v -= Y.dot((Xb.array() + model_param.intercept).matrix());

    for (int i = 0; i < n; i++)
      if (P[i] > 1e-8) v -= (log(P[i]) - model_param.intercept - Xb[i]);

    return (v / n);
  }

  PoissonObjective::PoissonObjective(const double *xmat, const double *y, int n,
                                     int d, int p, double step_size0, bool include_intercept)
    : GLMObjective(xmat, y, n, d, p, step_size0, include_intercept) {

    model_param.intercept = 0.0;
    update_auxiliary();
    for (int i = 0; i < d; i++) update_gradient(i);


    deviance = fabs(eval());
  }

  void PoissonObjective::update_auxiliary() {
    P = Xb.array() + model_param.intercept;
    P = P.array().exp();
    R = Y - P;
    W = P;
    sum_w = W.sum();
    wxx_gen++;
  }

  double PoissonObjective::eval() {
    double v = 0.0;
    for (int i = 0; i < n; i++)
      v = v + P[i] - Y[i] * (model_param.intercept + Xb[i]);
    return (v / n);
  }

}
