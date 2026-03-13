#ifndef SAM_ACTGD_H
#define SAM_ACTGD_H

#include "../objective/objective.h"
#include "solver_params.h"
#include <vector>

namespace SAM {
  class GroupActGDSolver {
  private:
    SolverParams m_param;
    ObjFunction *m_obj;
    std::vector<int> itercnt_path;

  public:
    std::vector<ModelParam> solution_path;
    GroupActGDSolver(ObjFunction *obj, SolverParams param);
    void solve(double *sse, int *df);

    ~GroupActGDSolver() {
      delete m_obj;
      m_obj = nullptr;
    }
  };
}  // namespace SAM

#endif
