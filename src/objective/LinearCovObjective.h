#ifndef LINEARCOV_OBJECTIVE_HPP
#define LINEARCOV_OBJECTIVE_HPP

#include "objective.h"

namespace SAM {

class LinearCovObjective : public ObjFunction {
private:
  MatrixXd C;      // (d*p) x (d*p): full Gram matrix X_all^T X_all / n
  VectorXd Xy;     // d*p: stacked X_j^T y / n
  VectorXd Xm;     // d*p: stacked column means X_j^T 1 / n
  double Ymean;    // mean of y
  double yy;       // ||y||^2 / n

public:
  LinearCovObjective(const double *xmat, const double *y, int n, int d, int p,
                     bool include_intercept);
  VectorXd coordinate_descent(RegFunction *regfunc, int idx);

  void intercept_update();
  void update_auxiliary();
  void update_gradient(int idx);

  double get_local_change(const VectorXd &old, int idx);
  double get_local_change_intercept(double old);

  double eval();
  double get_r2();
};

}  // namespace SAM

#endif
