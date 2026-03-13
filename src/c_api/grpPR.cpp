#include <Rcpp.h>
#include <RcppEigen.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include "../utils.h"
#include "../solver/actnewton.h"
#include "../objective/GLMObjective.h"
#include "c_api_utils.h"

using namespace SAM;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

extern "C" void grpPR(double *A, double* yy, double *lambda, int *nnlambda, double *LL0, int *nn, int *dd, int *pp, double *xx, double *aa0, int *mmax_ite, double *tthol, char** regfunc, double *aalpha, double *z, int *df, double *func_norm, int *ddfmax, int *vverbose, double *ddev_ratio_thr, double *ddev_change_thr) {

  double thol = *tthol, L0 = *LL0;
  int nlambda = *nnlambda, n = *nn, d = *dd, p = *pp, max_ite = *mmax_ite;

  VectorXd y(n);

  for (int i = 0; i < n; i++) {
    y(i) = yy[i];
  }
  for (int i = 0; i < nlambda; i++) {
    lambda[i] /= n;
  }

  SolverParams *param = new SolverParams();
  make_solver_params(param, lambda, nlambda, *regfunc, thol, max_ite, *ddfmax, *vverbose);
  param->dev_ratio_thr = *ddev_ratio_thr;
  param->dev_change_thr = *ddev_change_thr;

  ObjFunction *obj = new PoissonObjective(A, yy, n, d, p, L0, param->include_intercept);

  ActNewtonSolver solver(obj, *param);

  vector<double> sse(nlambda, 0.0);
  solver.solve(sse.data(), df);

  if (solver.solution_path.size() != (unsigned int)nlambda) {
    SAM_PRINTF("grpPR: solution_path size mismatch (%d vs %d)\n",
               (int)solver.solution_path.size(), nlambda);
    delete param;
    return;
  }
  for (int i = 0; i < nlambda; i++) {
    ModelParam &model = solver.solution_path[i];
    for (int j = 0; j < d; j++) {
      func_norm[i*d + j] = calc_norm(model.beta[j]);
      for (int k = 0; k < p; k++) {
        xx[i*(d*p+1) + j*p + k] = model.beta[j](k);
      }
    }
    xx[i*(d*p+1) + d*p] = model.intercept;
  }

  delete param;

}
