#ifndef SAM_SOLVER_PARAMS_H
#define SAM_SOLVER_PARAMS_H

#include <vector>

namespace SAM {
  enum RegType { L1, SCAD, MCP };

  // training parameters
  class SolverParams {
  public:
    /*! number of regularization parameters */
    unsigned num_lambda;

    /*! the last paramter on the regularization path */
    double target_lambda;

    /*! type of regularization terms */
    RegType reg_type;

    /*! gamma param for SCAD and MCP regularization */
    double gamma;

    /*！ rounds of relaxation when solving SCAD and MCP penalty */
    unsigned num_relaxation_round;

    /*! precision of optimization */
    double prec;

    /*! max number of iteration for innner loop */
    int max_iter;

    /*! whether or not to add intercept term */
    bool include_intercept;

    /*! maximum number of non-zero groups; -1 = unlimited */
    int dfmax;

    /*! print iteration info */
    bool verbose;

    /*! stop if deviance ratio (1 - cur_dev/null_dev) exceeds this; -1.0 = disabled */
    double dev_ratio_thr;

    /*! stop if relative deviance change < this over last min_lambda_count steps; -1.0 = disabled */
    double dev_change_thr;

    /*! minimum lambdas before checking dev_change_thr */
    int min_lambda_count;

    std::vector<double> lambdas;

    SolverParams();

    void set_lambdas(const double *lambda_path, int n);

    std::vector<double> get_lambda_path() const;
  };

}  // namespace SAM

#endif
