#ifndef SAM_C_API_UTILS_H
#define SAM_C_API_UTILS_H
#include "../solver/solver_params.h"
#include "../utils.h"
#include <cstring>

namespace SAM {

inline void make_solver_params(
    SolverParams *param,
    const double *lambda, int nlambda,
    const char *regfunc_str,
    double prec, int max_iter,
    int dfmax, int verbose)
{
    param->set_lambdas(lambda, nlambda);
    param->gamma = 3;
    if (std::strcmp(regfunc_str, "MCP") == 0)
        param->reg_type = MCP;
    else if (std::strcmp(regfunc_str, "SCAD") == 0)
        param->reg_type = SCAD;
    else
        param->reg_type = L1;
    param->include_intercept = true;
    param->prec = prec;
    param->max_iter = max_iter;
    param->num_relaxation_round = 10;
    param->dfmax = dfmax;
    param->verbose = (verbose != 0);
}

}  // namespace SAM
#endif
