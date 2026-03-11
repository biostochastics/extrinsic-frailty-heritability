# Heritable extrinsic susceptibility can inflate inferred intrinsic lifespan heritability

Simulation pipeline and manuscript demonstrating that **omitted familial extrinsic frailty** biases calibrate-then-extrapolate intrinsic heritability estimates upward by ~9-10 percentage points.

## Summary

[Shenhar et al. (2026, *Science*)](https://doi.org/10.1126/science.adz1187) estimate that intrinsic human lifespan heritability is ~50% by calibrating parametric mortality models to Scandinavian twin data and extrapolating to a counterfactual world without extrinsic mortality. [Hamilton (2026, *medRxiv*)](https://doi.org/10.1101/2026.01.15.26300001) argues that heritable infection mortality creates collider bias through concordant-survivor conditioning. However, Shenhar's procedure does not remove individuals --- it parametrically sets extrinsic hazard to zero.

We identify a distinct bias mechanism specific to Shenhar's calibrate-then-extrapolate approach: if susceptibility to extrinsic mortality is familial (including heritable) and this heterogeneity is not modeled, the single-parameter calibration absorbs extrinsic genetic variance into the intrinsic parameter --- omitted-variable bias at the estimation step, before any extrapolation occurs.

### Headline result (Gompertz-Makeham, 50-seed Monte Carlo)

| Condition | Bias (pp) | 95% CI |
|---|---|---|
| **Misspecified** (familial extrinsic omitted) | **+9.2** | 8.7 to 9.7 |
| Correctly specified (baseline) | -0.4 | -0.9 to 0.1 |

The fitted intrinsic frailty dispersion is inflated by ~21%, corresponding to ~9-10 pp of upward bias in inferred intrinsic h². The correctly specified baseline is indistinguishable from zero.

### Single-seed decomposition (illustrative)

| Condition | sigma\_theta,fit | h² | Bias (pp) |
|---|---|---|---|
| Oracle (true intrinsic) | 1.259 | 0.510 | --- |
| Correctly specified | 1.220 | 0.472 | -3.7 |
| **Misspecified** (sigma\_gamma=0.40, rho=0.4) | **1.537** | **0.586** | **+7.6** |
| Non-heritable extrinsic | 1.259 | 0.478 | -3.1 |
| Vanishing m\_ex | 1.267 | 0.500 | -0.9 |
| Irrelevant trait | 1.264 | 0.502 | -0.7 |
| Two-component recovery | 1.219 | 0.476 | -3.3 |
| Hamilton conditioning | N/A | 0.427 | -8.2 |

The net +7.6 pp decomposes into +11.4 pp inflation (omitted-variable absorption) partially offset by -3.7 pp attenuation (extrinsic noise).

### Cross-model replication (20 seeds per model)

| Model | Misspecified bias | 95% CI lower bound |
|---|---|---|
| Gompertz-Makeham | +9.2 +/- 0.3 pp | > 8.2 pp |
| Makeham-Gamma-Gompertz | +10.4 +/- 0.5 pp | > 8.2 pp |
| Saturating-Removal | +9.2 +/- 0.5 pp | > 8.2 pp |

All three mortality-model families show the same phenomenon with all 95% CIs above +8.2 pp.

### Anchored sensitivity regime

Parameters anchored to independent immunogenetic evidence: sigma\_gamma from twin heritability of infection mortality ([Obel et al. 2010](https://doi.org/10.1093/aje/kwq037)) and rho from LDSC genetic correlations between infection severity and longevity ([Qiu et al. 2025](https://doi.org/10.1186/s12967-024-05932-y)).

| Region | Bias range (pp) | Mean |
|---|---|---|
| Anchored regime (sigma\_gamma in [0.30, 0.65], rho in [0.20, 0.50]) | +2.8 to +18.9 | +9 |

### Robustness

- **Functional forms**: Log-normal (+8.1 pp), additive (+8.7 pp), gamma (+7.6 pp) extrinsic frailty constructions (10 seeds each)
- **Sign reversal**: Negative intrinsic-extrinsic genetic correlation reverses the bias direction (29 rho values x 20 seeds, three model families)
- **Bivariate diagnostic**: Misspecification produces structured age-specific misfit in twin dependence (IAE ratio 4.7x), removed by correct specification
- **Multi-target calibration**: Richer targets reduce but do not eliminate inflation (+10.0 pp -> +6.4 pp)

## Repository structure

```
├── writeup-biorxiv-v3.pdf     # Compiled manuscript
├── writeup-biorxiv-v3.typ     # Typst source (dynamic data binding from scalars.json)
├── proof-draft.md             # Formal proof (Appendix B source)
├── _targets.R                 # Pipeline definition (194 targets, static branching)
├── R/
│   ├── params.R               # Single source of truth for all parameters
│   ├── sim_core.R             # Lambert W GM solver, hybrid MGG solver, SR Euler-Maruyama
│   ├── sim_twins.R            # Twin parameter generation, sim_twin_h2()
│   ├── calibrate.R            # Stochastic bisection (CRN), oracle, joint, multi-target, UQ
│   ├── utils.R                # Sweep cells, Hamilton arm, model analysis, control sweeps
│   ├── anchors.R              # sigma_gamma bridge, tetrachoric, m_ex split, negative rho
│   ├── diagnostics.R          # Variance decomposition, joint r_MZ/r_DZ, alpha/beta estimation
│   ├── replication.R          # 50-seed Monte Carlo replication
│   ├── plots.R                # All figure-generating functions
│   ├── commentary_figure.R    # Science commentary figure (Fig. 1, 5 panels)
│   └── export.R               # CSV/JSON export
├── src/
│   └── sim_lifespan_sr.cpp    # Rcpp/C++ SR solver (Xoroshiro128+ RNG, ~4.6x speedup)
├── tables/                    # Pipeline outputs (CSV/JSON)
│   ├── scalars.json           # All scalar results (dynamic data binding for Typst)
│   ├── summary.csv            # Main summary table
│   ├── models.csv             # Cross-model validation
│   ├── sweep.csv              # Full 81-cell sensitivity grid
│   ├── anchored.csv           # 56-cell anchored grid
│   ├── controls.csv           # Control experiments
│   └── model_controls.csv     # MGG/SR control sweeps
└── figures/                   # All generated figures (PNG, 29 files)
```

## Reproducing

### Requirements

**R** (>= 4.4) with packages:
```r
install.packages(c("targets", "tarchetypes", "crew",
                    "MASS", "ggplot2", "patchwork", "lamW",
                    "ggpubr", "ggthemes", "viridis", "mvtnorm",
                    "scales", "Rcpp"))
```

### Run the pipeline

```bash
git clone https://github.com/biostochastics/extrinsic-frailty-heritability.git
cd extrinsic-frailty-heritability

# Full pipeline (~20 hours wall-clock on Apple M4 Pro, 10 parallel workers via crew)
Rscript -e 'targets::tar_make()'

# Inspect pipeline DAG
Rscript -e 'targets::tar_visnetwork()'

# Read a specific target
Rscript -e 'targets::tar_read(sweep_results)'
```

The pipeline uses `{crew}` for parallel execution with 10 local workers. It defines 194 targets organized in a DAG with `tar_map` static branching for the 81-cell sensitivity sweep and 56-cell anchored sweep. A master seed (66829) propagates deterministic sub-seeds to every target for full reproducibility.

### Pipeline stages

```
Stage 0:  Calibration          -> sigma_theta_true, oracle
Stage 1:  Main arms            -> Arms 1-3, Controls A/B, Oracle fix
Stage 2:  Sensitivity sweep    -> 81 cells (9x9 tar_map static branching)
Stage 3:  Anchored sweep       -> 56 cells (7x8 tar_map static branching)
Stage 4:  Variance decomp      -> alpha, beta coefficients (R^2 = 0.988)
Stage 5:  Model validation     -> MGG + SR families
Stage 6:  Additional controls  -> Dose-response, pleiotropy isolation, irrelevant trait
Stage 7:  Diagnostics          -> Joint r_MZ/r_DZ, sigma_gamma bridge, negative rho
Stage 8:  Calibration UQ       -> Bootstrap uncertainty
Stage 8b: Revision analyses    -> Functional forms (RC1), extended sigma_gamma (RC2),
                                  MC uncertainty (RC3/RC3b), heritable fraction (RC4),
                                  bivariate check (B5), MZ concordance (B6/B6b),
                                  bridge uncertainty (B7)
Stage 9:  Tables               -> CSV/JSON export
Stage 10: Figures              -> All plots (29 figures)
```

## Technical notes

- **GM solver**: Closed-form via Lambert W function (`lamW`), with bisection fallback.
- **MGG parameterization**: `a = a_0^q` (Strehler-Mildvan compensation). Shenhar's supplement states `a = a_0^(1/q)`, which does not produce hazard curve crossing and yields inconsistent twin correlations.
- **SR solver**: Euler-Maruyama at weekly resolution (dt = 1/52), compiled to C++ via Rcpp with Xoroshiro128+ RNG (~4.6x speedup over vectorized R).
- **Calibration**: 40-iteration stochastic bisection with common random numbers (CRN) for variance reduction.
- **50,000 twin pairs** per zygosity per simulation (GM/MGG); 25,000 for SR.

## Manuscript

The compiled manuscript (`kornilov_manuscript_extended.pdf`) includes the full technical report with two appendices:
- **Appendix A** (A.1-A.12): Computational implementation details
- **Appendix B** (B.1-B.9): Formal proof of calibration inflation (delta-method analysis)

Scalar results in `tables/scalars.json` are dynamically bound into the Typst source at compile time.

## Author

Sergey A. Kornilov, Biostochastics, Seattle, WA

## Acknowledgements

GPT-5.2-Pro (OpenAI) and Gemini-3-Pro (Google) were used during development to review the simulation design, verify that the calibration procedure parallels Shenhar et al.'s methodology, and provide structural feedback. MiniMax-M2.5 (MiniMax), GLM-5 (Zhipu AI), GPT-5.2-Pro, and Gemini-3-Pro were used to review and assist in deriving the formal proofs. Claude Opus 4.6 (Anthropic) was used for code generation, drafting assistance, and proof typesetting. All code, analytical results, and final text were reviewed and revised by the author.

## License

MIT
