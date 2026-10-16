# Sieve-free ReLU Bayesian priors with near-minimax posterior contraction

R code and simulation results for the [paper](paper/main_arxiv.pdf) by Hien Duy Nguyen, Jacob Westerhout and Julyan Arbel.

The study compares Gaussian, uniform and Student parameter priors for two regression targets, with likelihood power `eta = 0.9`, sample sizes from 64 to 4096 and 50 replications per case.

## Requirements

The reference simulations use R 4.4.0. Fitting requires a C++ compiler configured for R and the following packages:

```r
install.packages(c("Rcpp", "RcppArmadillo", "coda"))
```

Recreating the figures also requires Python 3 with `pypdf`. Poppler's `pdftoppm` is optional for PNG output.

## Reproduce and verify

Run from the repository root:

```sh
Rscript Replot.R unused --check-only
Rscript Refit.R unused --check-only
Rscript Replot.R reconstructed
```

The first two commands check the supplied results and simulation settings. The third recreates the figures and tables from the included numerical results, saving them in `reconstructed/`.

To refit the reported models:

```sh
Rscript Refit.R refitted --workers=4
```

Refitting the complete study is computationally expensive. Add `--rows=1` to refit one case, or `--all-stages` to follow the computation stages from the beginning. Choose a new output directory for each run. The seeds, selected settings and diagnostic requirements are included in `records/` and `protocol.json`.

## Source files

| File | Purpose |
| --- | --- |
| `Refit.R` | Refit the simulation models |
| `Replot.R` | Check results and recreate figures and tables |
| `source/` | R functions and C++ sampling code |
| `records/` | Replication results, seeds and fitted settings |
| `protocol.json` | Simulation settings and diagnostic requirements |
| `paper/` | Article source and PDF |

## Reference outputs

The included results contain 2,100 fitted models and 70 mean squared-error comparisons, each based on 50 replications.

| Manuscript output | File |
| --- | --- |
| Figure 1 | `figures/uniform_lipschitz_credible_intervals.pdf` |
| Figure 2 | `figures/adaptive_squared_error.pdf` |
| Figure 3 | `figures/two_lobes_credible_intervals.pdf` |
| Tables 3-7 | `figures/*_table.tex` |

## Contact

For questions about the code, please use the [repository issues](https://github.com/hiendn/bayesian-relu-minimax/issues).
