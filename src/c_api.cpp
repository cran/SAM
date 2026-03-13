#include "c_api/grplasso.cpp"
#include "c_api/grpSVM.cpp"
#include "c_api/grpLR.cpp"
#include "c_api/grpPR.cpp"
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


static R_NativePrimitiveArgType grplasso_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, STRSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, INTSXP
};
static R_NativePrimitiveArgType grpSVM_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP
};
static R_NativePrimitiveArgType grpLR_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, STRSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP
};
static R_NativePrimitiveArgType grpPR_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, STRSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP
};
static const R_CMethodDef c_Methods[] = {
  {"grplasso", (DL_FUNC) &grplasso, 21, grplasso_t},
  {"grpSVM", (DL_FUNC) &grpSVM, 18, grpSVM_t},
  {"grpLR", (DL_FUNC) &grpLR, 21, grpLR_t},
  {"grpPR", (DL_FUNC) &grpPR, 21, grpPR_t},
  {NULL, NULL, 0, NULL}
};
extern "C" void R_init_SAM(DllInfo *info)
{
  R_registerRoutines(info, c_Methods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
