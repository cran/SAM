#ifndef UTILS_HPP
#define UTILS_HPP

#include <Rcpp.h>
#include <RcppEigen.h>
using Eigen::VectorXd;
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::plugins(openmp)]

// Portable printf: Rprintf for R, printf for standalone (Python ctypes)
#ifdef SAM_STANDALONE
  #include <cstdio>
  #define SAM_PRINTF std::printf
  #define SAM_FFLUSH() std::fflush(stdout)
#else
  #define SAM_PRINTF Rprintf
  #define SAM_FFLUSH() ((void)0)
#endif

namespace SAM {
  extern double calc_norm(const VectorXd &x);
  extern double sqr(double x);
}

#endif
