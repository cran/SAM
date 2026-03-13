#include <Rcpp.h>
#include <RcppEigen.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include "../utils.h"
#include "../solver/actnewton.h"
#include "../solver/actgd.h"
#include "../objective/LinearObjective.h"
#include "../objective/LinearCovObjective.h"
#include "c_api_utils.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]


using namespace SAM;

extern "C" void grplasso(double *yy, double *XX, double *lambda, int *nnlambda, int *nn, int *dd, int *pp, double *ww, int *mmax_ite, double *tthol, char** regfunc, int *iinput, int *df, double *sse, double *func_norm, int *ddfmax, int *vverbose, double *ddev_ratio_thr, double *ddev_change_thr, int *ssolver_type, int *ttype_gaussian)
{
  int n,d,p,max_ite,nlambda;
  int input;

  double thol;
  double lambda_max;

  nlambda = *nnlambda;
  n = *nn;
  d = *dd;
  p = *pp;
  max_ite = *mmax_ite;
  thol = *tthol;
  input = *iinput;

  lambda_max = 0;

  vector<MatrixXd> V(d);
  VectorXd y(n);
  const double svd_tol = 1e-12;

  for (int i = 0; i < n; i++)
    y(i) = yy[i];


  for (int i = 0; i < d; i++){

    MatrixXd X(n, p);
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < p; k++) {
        X(j,k) = XX[i*n*p + k*n + j];
      }
    }

    //std::cout << "X:" << std::endl << X << std::endl;

    Eigen::JacobiSVD<MatrixXd> svd(X, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd S = svd.singularValues();

    X = svd.matrixU();

    V[i] = svd.matrixV();
    for (int j = 0; j < p; j++)
      for (int k = 0; k < p; k++)
        V[i](j,k) /= std::max(S[k], svd_tol);

    for (int j = 0; j < n; j++) {
      for (int k = 0; k < p; k++) {
        XX[i*n*p + k*n + j] = X(j, k);
      }
    }

    if (input == 0) {
      lambda_max = std::max(lambda_max, calc_norm(X.transpose()*y)/n);
    }
  }

  if (input == 0) {
    for (int i = 0; i < nlambda; i++)
      lambda[i] = lambda[i] * lambda_max;
  }

  SolverParams *param = new SolverParams();
  make_solver_params(param, lambda, nlambda, *regfunc, thol, max_ite, *ddfmax, *vverbose);
  param->dev_ratio_thr = *ddev_ratio_thr;
  param->dev_change_thr = *ddev_change_thr;

  // type_gaussian: 0 = naive, 1 = covariance (fallback to naive if d*p >= n),
  //                2 = auto (covariance if d*p < n, else naive)
  bool use_cov = (*ttype_gaussian == 1 || *ttype_gaussian == 2) && d * p < n;
  ObjFunction *obj;
  if (use_cov) {
    obj = new LinearCovObjective(XX, yy, n, d, p, param->include_intercept);
  } else {
    obj = new LinearObjective(XX, yy, n, d, p, param->include_intercept);
  }

  // solver dispatch: 1 = ActGD (L1 only), 0 = ActNewton (default)
  if (*ssolver_type == 1 && param->reg_type == L1) {
    GroupActGDSolver solver(obj, *param);
    solver.solve(sse, df);
    if (solver.solution_path.size() != (unsigned int)nlambda) {
      SAM_PRINTF("grplasso: solution_path size mismatch (%d vs %d)\n",
                 (int)solver.solution_path.size(), nlambda);
      delete param;
      return;
    }
    // extract results
    for (int i = 0; i < nlambda; i++) {
      ModelParam &model = solver.solution_path[i];
      for (int j = 0; j < d; j++) {
        func_norm[i*d+j] = calc_norm(model.beta[j]);
        VectorXd res = V[j] * model.beta[j];
        for (int k = 0; k < p; k++) {
          ww[i*d*p + j*p + k] = res(k);
        }
      }
    }
  } else {
    ActNewtonSolver solver(obj, *param);
    solver.solve(sse, df);
    if (solver.solution_path.size() != (unsigned int)nlambda) {
      SAM_PRINTF("grplasso: solution_path size mismatch (%d vs %d)\n",
                 (int)solver.solution_path.size(), nlambda);
      delete param;
      return;
    }
    for (int i = 0; i < nlambda; i++) {
      ModelParam &model = solver.solution_path[i];
      for (int j = 0; j < d; j++) {
        func_norm[i*d+j] = calc_norm(model.beta[j]);
        VectorXd res = V[j] * model.beta[j];
        for (int k = 0; k < p; k++) {
          ww[i*d*p + j*p + k] = res(k);
        }
      }
    }
  }

  delete param;
}
