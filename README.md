# SAM

Sparse Additive Modelling (SAM) with:

- R package implementation (`R/`, `src/`, `man/`)
- Python wrapper package (`python-package/`) that reuses the current C++ core

## R package

Build/check:

```bash
R CMD build .
R CMD check --as-cran SAM_1.2.tar.gz
```

## Python package (wrapper)

Location: `python-package/`

APIs implemented:

- `samLL`, `samHL`, `samEL`, `samQL`
- `predict_samLL`, `predict_samHL`, `predict_samEL`, `predict_samQL`
- path plotting and summary helpers

The wrapper calls native symbols from `src/SAM.so`:

- `grpLR`, `grpSVM`, `grpPR`, `grplasso`

Quick local run:

```bash
bash python-package/scripts/build_native.sh
PYTHONPATH=python-package python3 python-package/examples/run_smoke.py
```

## Documentation

- R help files: `man/*.Rd`
- Python docs: `python-package/docs/` + `python-package/mkdocs.yml`

To build Python docs locally:

```bash
cd python-package
mkdocs build --strict
```
