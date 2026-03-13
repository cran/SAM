#include <Rcpp.h>
#include <RcppEigen.h>
#include "actgd.h"
#include "../objective/objective.h"
#include "solver_params.h"
#include "../utils.h"
#include <algorithm>
#include <vector>

//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

namespace SAM {

GroupActGDSolver::GroupActGDSolver(ObjFunction *obj, SolverParams param)
    : m_param(param), m_obj(obj) {
  itercnt_path.clear();
  solution_path.clear();
}

void GroupActGDSolver::solve(double *sse, int *df) {
  int d = m_obj->get_dim();

  const std::vector<double> &lambdas = m_param.get_lambda_path();
  itercnt_path.resize(lambdas.size(), 0);

  double dev_thr = m_obj->get_deviance() * m_param.prec;

  // strong_set[j] == 1 means group j passed strong rule screening
  std::vector<int> strong_set(d, 0);
  std::vector<int> strong_set_master(d, 0);
  std::vector<double> grad_norms(d, 0.0);
  std::vector<double> grad_norms_master(d, 0.0);

  // L1 regularization only
  RegFunction *regfunc = new RegL1();

  // compute initial gradient norms
  m_obj->update_auxiliary();
  for (int j = 0; j < d; j++) {
    grad_norms[j] = calc_norm(m_obj->get_grad(j));
  }

  // master model for warm-start (mirrors ActNewton's model_master)
  ModelParam model_master = m_obj->get_model_param();
  VectorXd Xb_master = m_obj->get_model_Xb();
  for (int j = 0; j < d; j++) grad_norms_master[j] = grad_norms[j];

  std::vector<double> dev_history(lambdas.size(), 0.0);
  std::vector<VectorXd> old_coef(d);

  for (size_t i = 0; i < lambdas.size(); i++) {
    double lambda_i = lambdas[i];

    // restore from master (same as ActNewton)
    m_obj->set_model_param(model_master);
    m_obj->set_model_Xb(Xb_master);
    for (int j = 0; j < d; j++) {
      grad_norms[j] = grad_norms_master[j];
      strong_set[j] = strong_set_master[j];
    }

    // strong rule screening
    double strong_thr;
    if (i > 0)
      strong_thr = 2.0 * lambda_i - lambdas[i - 1];
    else
      strong_thr = 2.0 * lambda_i;

    for (int j = 0; j < d; j++) {
      if (grad_norms[j] > strong_thr)
        strong_set[j] = 1;
    }

    m_obj->update_auxiliary();

    // Level 1 loop: active set update (mirrors ActNewton's loop level 1)
    int level1_iter = 0;
    std::vector<int> active_idx;

    while (level1_iter < m_param.max_iter) {
      level1_iter++;

      double old_intcpt = m_obj->get_intercept();
      for (int j = 0; j < d; j++) old_coef[j] = m_obj->get_model_coef(j);

      // one pass of CD on strong_set → identify active_idx (non-zero groups)
      active_idx.clear();
      for (int j = 0; j < d; j++) {
        if (!strong_set[j]) continue;
        regfunc->set_param(lambda_i, 0.0);
        VectorXd updated = m_obj->coordinate_descent(regfunc, j);
        if (calc_norm(updated) > 0)
          active_idx.push_back(j);
      }

      // Level 2 loop: iterate CD on active_idx until convergence
      int level2_iter = 0;
      while (level2_iter < m_param.max_iter) {
        level2_iter++;
        bool level2_converged = true;

        for (int k = 0; k < (int)active_idx.size(); k++) {
          int idx = active_idx[k];
          VectorXd old_beta = m_obj->get_model_coef(idx);
          regfunc->set_param(lambda_i, 0.0);
          m_obj->coordinate_descent(regfunc, idx);
          if (m_obj->get_local_change(old_beta, idx) > dev_thr)
            level2_converged = false;
        }

        if (m_param.include_intercept) {
          double oi = m_obj->get_intercept();
          m_obj->intercept_update();
          if (m_obj->get_local_change_intercept(oi) > dev_thr)
            level2_converged = false;
        }

        if (level2_converged) break;
      }

      itercnt_path[i] += level2_iter;

      // check Level 1 convergence
      bool level1_converged = true;
      for (int k = 0; k < (int)active_idx.size(); k++) {
        int idx = active_idx[k];
        if (m_obj->get_local_change(old_coef[idx], idx) > dev_thr)
          level1_converged = false;
      }
      if (m_param.include_intercept &&
          m_obj->get_local_change_intercept(old_intcpt) > dev_thr)
        level1_converged = false;

      m_obj->update_auxiliary();

      if (!level1_converged) continue;

      // Active set converged — check KKT on all groups outside strong_set
      bool new_violations = false;
      for (int j = 0; j < d; j++) {
        if (strong_set[j] == 0) {
          m_obj->update_gradient(j);
          double gnorm = calc_norm(m_obj->get_grad(j));
          grad_norms[j] = gnorm;
          if (gnorm > lambda_i) {
            strong_set[j] = 1;
            new_violations = true;
          }
        }
      }

      // No KKT violations — done for this lambda
      if (!new_violations) break;
      // else: violations found, continue with expanded strong set
    }

    // save master state (mirrors ActNewton's model_master save after level_0 == 1)
    {
      const ModelParam &ref = m_obj->get_model_param_ref();
      const VectorXd &Xb_ref = m_obj->get_model_Xb_ref();
      model_master.intercept = ref.intercept;
      for (int j = 0; j < d; j++) {
        model_master.beta[j] = ref.beta[j];
        grad_norms_master[j] = grad_norms[j];
        strong_set_master[j] = strong_set[j];
      }
      Xb_master = Xb_ref;
    }

    solution_path.push_back(m_obj->get_model_param());

    // df = number of groups in active_idx (matches ActNewton)
    df[i] = (int)active_idx.size();
    sse[i] = m_obj->get_r2();

    double null_dev = m_obj->get_null_deviance();
    double cur_dev = m_obj->get_current_deviance();
    dev_history[i] = cur_dev;

    if (m_param.verbose) {
      SAM_PRINTF("lambda %d/%d: df=%d, iter=%d\n",
              (int)(i + 1), (int)lambdas.size(), df[i], itercnt_path[i]);
      SAM_FFLUSH();
    }

    // early stopping: dfmax
    if (m_param.dfmax > 0 && df[i] >= m_param.dfmax) {
      for (size_t j = i + 1; j < lambdas.size(); j++) {
        solution_path.push_back(m_obj->get_model_param());
        df[j] = df[i]; sse[j] = sse[i]; itercnt_path[j] = 0;
      }
      break;
    }

    // early stopping: deviance ratio
    if (m_param.dev_ratio_thr > 0 && df[i] > 0 && null_dev > 0) {
      double dev_ratio = 1.0 - cur_dev / null_dev;
      if (dev_ratio > m_param.dev_ratio_thr) {
        if (m_param.verbose) {
          SAM_PRINTF("  early stop: dev_ratio=%.6f > %.6f\n",
                  dev_ratio, m_param.dev_ratio_thr);
          SAM_FFLUSH();
        }
        for (size_t j = i + 1; j < lambdas.size(); j++) {
          solution_path.push_back(m_obj->get_model_param());
          df[j] = df[i]; sse[j] = sse[i]; itercnt_path[j] = 0;
        }
        break;
      }
    }

    // early stopping: relative deviance change
    if (m_param.dev_change_thr > 0 && (int)i >= m_param.min_lambda_count && df[i] > 0) {
      double prev_dev = dev_history[i - m_param.min_lambda_count];
      if (fabs(prev_dev) > 1e-10) {
        double rel_change = fabs(cur_dev - prev_dev) / fabs(prev_dev);
        if (rel_change < m_param.dev_change_thr) {
          if (m_param.verbose) {
            SAM_PRINTF("  early stop: dev_change=%.2e < %.2e\n",
                    rel_change, m_param.dev_change_thr);
            SAM_FFLUSH();
          }
          for (size_t j = i + 1; j < lambdas.size(); j++) {
            solution_path.push_back(m_obj->get_model_param());
            df[j] = df[i]; sse[j] = sse[i]; itercnt_path[j] = 0;
          }
          break;
        }
      }
    }
  }

  delete regfunc;
}

}  // namespace SAM
