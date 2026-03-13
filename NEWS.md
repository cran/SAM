# SAM Changelog

## v1.3

### Bug Fixes
- Fixed KKT validation: inactive groups are now checked after every Level 1
  iteration, not only when not converged.
- Fixed `grpSVM` `ya_idx` active-set logic that only checked the first element
  of each group instead of all `p` coefficients.
- Replaced `assert()` with runtime error checks in `grplasso`, `grpLR`, and
  `grpPR` (asserts are compiled out under `-DNDEBUG`).
- Added NULL checks for all `malloc` calls in `grpSVM`.

### New Features
- `dfmax` parameter for early stopping when the number of non-zero groups
  exceeds a threshold (all four families).
- `verbose` parameter to print per-lambda iteration info.
- Deviance-ratio and deviance-change early stopping (`dev.ratio.thr`,
  `dev.change.thr`) for samQL, samLL, and samEL.
- Unified `sam()` entry point dispatching by `family` argument (R and Python).
- `coef()` S3 methods for all four model classes.
- `solver = "actgd"` option for samQL (group active gradient descent with
  strong-rule screening, L1 only).
- `type.gaussian` option for samQL (`"naive"` / `"covariance"` / `"auto"`).
- `runtime` field on fitted objects (R: `$runtime`, Python: `.runtime_seconds`).

### Python Wrapper
- New `pysam-sam` Python package calling the same C++ core via `ctypes`.
- Pure-Python spline basis (no R dependency at runtime).
- Full API parity: `samLL`, `samHL`, `samEL`, `samQL` + predict functions.

### Performance
- Cached `wXX` products and no-op `update_auxiliary` for linear objective.
- C API deduplication via `c_api_utils.h`.

### Code Quality
- Removed redundant `fit$nkots` field from R wrappers.
- Removed unused `fit$XX` (full design matrix) from samQL output.
- Removed unused `V(d)` declarations in `grpLR` and `grpPR`.
- Extracted shared `.sam_print_path()` / `.sam_plot_path()` helpers.
- Removed unused `scipy` dependency from Python package.
- Fixed legacy `PICASSO_SOLVER_PARAMS_H` header guard.
- Added root `.gitignore`.
- Added test coverage for samHL and samEL.

## v1.2

- Initial CRAN release with samQL, samLL, samEL, samHL.
