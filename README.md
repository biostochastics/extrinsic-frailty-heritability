# Omitted familial extrinsic risk inflates inferred intrinsic lifespan heritability

Simulation pipeline demonstrating that **omitted familial extrinsic frailty** biases calibrate-then-extrapolate intrinsic heritability estimates upward. Response to Shenhar et al. (2026, *Science*) and Hamilton (2026, *medRxiv*).

## Summary

[Shenhar et al. (2026, *Science*)](https://doi.org/10.1126/science.adz1187) estimate ~50% intrinsic human lifespan heritability by calibrating parametric mortality models to Scandinavian twin data and extrapolating to a counterfactual world without extrinsic mortality. We identify a bias mechanism specific to their calibrate-then-extrapolate approach: if susceptibility to extrinsic mortality is familial (heritable) and this heterogeneity is not modeled, the single-parameter calibration absorbs extrinsic genetic variance into the intrinsic frailty parameter — omitted-variable bias at the estimation step, before any extrapolation occurs.

### Two-layer vulnerability

The bias operates at two separable levels:

1. **Model-class problem:** Even full maximum likelihood estimation within the misspecified one-component model inflates sigma_theta by +3–5%, because the model has no parameter for the omitted extrinsic component.
2. **Estimator-choice problem:** Shenhar's moment-based calibration (matching r_MZ) amplifies the inflation ~7×, producing +21–42% inflation, because the omitted pathway loads directly into twin dependence and the calibrator has no marginal-survival counterweight.

### Primary estimand: sigma_theta inflation

| Model | sigma\_true | sigma\_fit | sigma inflation | Falconer bias (pp) |
|---|---|---|---|---|
| **Gompertz-Makeham** | 1.259 | 1.537 | **+22.1%** | +7.6 |
| MGG | 0.270 | 0.318 | +17.9% | +9.2 |
| SR | 3.570 | 4.073 | +14.1% | +5.7 |
| Recovery (GM) | 1.259 | 1.219 | -3.2% | -3.3 |

Monte Carlo uncertainty (20 seeds): sigma_theta inflation +22.1 +/- 0.3% SE (95% CI: 21.5–22.7%).

### ML vs moment calibration

| Condition | ML inflation | Moment inflation | Ratio |
|---|---|---|---|
| Correctly specified (sigma\_gamma=0) | +0.0% | ~0% | — |
| Default (sigma\_gamma=0.40, rho=0.4) | **+3.0%** | **+21.3%** | **~7×** |
| No pleiotropy (rho=0) | +0.3% | +8.8% | ~29× |
| Extreme (sigma\_gamma=0.65) | +4.8% | +41.8% | ~9× |

Under misspecification, ML converges to a different pseudo-true parameter than moment calibration (White 1982): ML minimizes KL divergence over the full bivariate density, moment calibration reproduces r_MZ. The omitted extrinsic pathway loads into twin dependence, so ML — partially regularized by marginal survival structure — absorbs far less of the omitted variance.

### Anchored sensitivity regime

Parameters anchored to independent immunogenetic evidence:

| Region | Bias range (pp) | Mean |
|---|---|---|
| Anchored (sigma\_gamma in [0.30, 0.65], rho in [0.20, 0.50]) | +2.8 to +18.9 | +9 |

### Evidence hierarchy

**Tier 1 — Core argument:** sigma_theta inflation (+22.1%), two-component decomposition and recovery (-3.2%), variance absorption diagnostic (R² = 0.988), independence from survival conditioning.

**Tier 2 — Structural fingerprints:** Bivariate survival surface misfit (ISE 48×), age-specific conditional correlation (peak Delta r = 0.048 at age 80), ACE decomposition (C ~ 0; bias loads onto A, not shared environment).

**Tier 3 — Scope, specificity, falsifiability:** Negative controls (3 necessary conditions), model generality (GM/MGG/SR), dose-response, anchored quantification, sign reversal under negative rho (confirmed in all 3 models), ML vs moment comparison.

## Repository structure

```
├── _targets.R                 # Pipeline definition (210 targets, static branching)
├── R/
│   ├── params.R               # Single source of truth for all parameters
│   ├── sim_core.R             # Lambert W GM solver, hybrid MGG solver, SR Euler-Maruyama
│   ├── sim_twins.R            # Twin parameter generation, sim_twin_h2(), sim_twin_lifespans()
│   ├── calibrate.R            # Stochastic bisection (CRN), oracle, joint, multi-target, UQ
│   ├── ml_estimation.R        # Full bivariate ML (GH quadrature, truncation/censoring)
│   ├── utils.R                # Sweep cells, Hamilton arm, model analysis, control sweeps
│   ├── anchors.R              # sigma_gamma bridge, tetrachoric, m_ex split, negative rho
│   ├── diagnostics.R          # Variance decomposition, joint r_MZ/r_DZ, alpha/beta estimation
│   ├── replication.R          # 50-seed Monte Carlo replication
│   ├── plots.R                # All figure-generating functions (33 plot functions)
│   ├── export.R               # CSV/JSON export
│   ├── variance_components.R  # sigma_theta inflation reporting, sweep decomposition, MC UQ
│   ├── bivariate_survival.R   # Joint S(t1,t2), ISE/sup/CvM/IAE metrics, concordance C(t)
│   ├── ace_decomp.R           # ACE/ADE twin decomposition via OpenMx
│   ├── age_dependence.R       # Conditional correlation, Kendall tau, cross-ratio by age
│   └── export_estimands.R     # Standalone CSV/JSON export for estimand analyses
├── src/
│   └── sim_lifespan_sr.cpp    # Rcpp/C++ SR solver (Xoroshiro128+ RNG, ~4.6x speedup)
├── tables/                    # Pipeline outputs (CSV/JSON)
│   ├── scalars.json           # All scalar results
│   ├── ml_scalars.json        # ML vs moment comparison results
│   ├── ml_vs_moment.csv       # Full ML experiment summary
│   └── ...                    # Summary, sweep, anchored, controls, models CSVs
├── figures/                   # All generated figures (43 PNG files)
│   ├── fig_ml_scatter.png     # ML vs moment scatter
│   ├── fig_ml_convergence.png # ML convergence strip plot
│   └── ...                    # Main results, robustness, diagnostics, estimand figures
├── lean/                      # Lean 4 formal proofs of Appendix B propositions
│   ├── CalibrationInflation.lean          # Root import (all modules)
│   ├── CalibrationInflation/
│   │   ├── Params.lean                    # Parameter structure + positivity lemmas
│   │   ├── Monotonicity.lean              # Prop B2: Falconer h² strict monotonicity
│   │   ├── Moments.lean                   # Moment definitions, feasibility condition
│   │   ├── Calibration.lean               # Calibration equation existence + expanded form
│   │   ├── Inflation.lean                 # Prop B1: variance inflation (ρ ≥ 0)
│   │   └── SignReversal.lean              # Cor B3: deflation under negative pleiotropy
│   ├── lakefile.toml                      # Lake build configuration (Mathlib dependency)
│   └── lean-toolchain                     # Lean 4 v4.29.0-rc8
└── webapp/                    # Interactive visualizations
```

## Reproducing

### Requirements

**R** (>= 4.4) with packages:
```r
install.packages(c("targets", "tarchetypes", "crew",
                    "MASS", "ggplot2", "patchwork", "lamW",
                    "ggpubr", "ggthemes", "viridis", "mvtnorm",
                    "scales", "Rcpp", "OpenMx", "statmod", "jsonlite"))
```

### Run the pipeline

```bash
git clone https://github.com/biostochastics/extrinsic-frailty-heritability.git
cd extrinsic-frailty-heritability

# Full pipeline (~20 hours wall-clock on Apple M4 Pro, 10 parallel workers via crew)
Rscript -e 'targets::tar_make()'

# ML comparison only (~4 minutes)
Rscript -e 'targets::tar_make(names = c("ml_experiment", "ml_figures", "ml_export"))'
```

The pipeline uses `{crew}` for parallel execution with 10 local workers. A master seed (66829) propagates deterministic sub-seeds to every target for full reproducibility.

### Pipeline stages

```
Stage 0:  Calibration          -> sigma_theta_true, oracle
Stage 1:  Main arms            -> Arms 1-3, Controls A/B, Oracle fix
Stage 2:  Sensitivity sweep    -> 81 cells (9x9 tar_map static branching)
Stage 3:  Anchored sweep       -> 56 cells (7x8 tar_map static branching)
Stage 4:  Variance decomp      -> alpha, beta coefficients (R² = 0.988)
Stage 5:  Model validation     -> MGG + SR families
Stage 6:  Additional controls  -> Dose-response, pleiotropy isolation, irrelevant trait
Stage 7:  Diagnostics          -> Joint r_MZ/r_DZ, sigma_gamma bridge, negative rho
Stage 8:  Calibration UQ       -> Bootstrap uncertainty
Stage 8b: Revision analyses    -> Functional forms, extended sigma_gamma, MC uncertainty,
                                  heritable fraction, bivariate check, MZ concordance
Stage 9:  Tables               -> CSV/JSON export
Stage 10: Figures              -> All plots
Stage 11: Alternative estimands -> Variance components, bivariate survival, ACE, age dependence
Stage 12: ML comparison        -> Full bivariate ML vs moment calibration (4 conditions × 3 sizes)
```

## Technical notes

- **GM solver**: Closed-form via Lambert W function (`lamW`), with bisection fallback.
- **MGG parameterization**: `a = a_0^q` (Strehler-Mildvan compensation).
- **SR solver**: Euler-Maruyama at weekly resolution (dt = 1/52), compiled to C++ via Rcpp with Xoroshiro128+ RNG (~4.6× speedup).
- **ML estimation**: Bivariate MZ log-likelihood with Gauss-Hermite quadrature (64 nodes), left-truncation correction (age-15 cutoff), right-censoring at t_max=120, mean-frailty constraint E[A]=const.
- **Calibration**: 40-iteration stochastic bisection with CRN for variance reduction.
- **50,000 twin pairs** per zygosity per simulation (GM/MGG); 25,000 for SR.

## Formal verification (Lean 4)

The `lean/` directory contains machine-checked proofs of the analytic results in Appendix B using [Lean 4](https://lean-lang.org/) with [Mathlib](https://leanprover-community.github.io/). The formalization covers:

- **Proposition B1** — Under non-negative pleiotropy (ρ ≥ 0) and feasibility, the calibrated intrinsic dispersion exceeds the true value: ŝ_θ > σ_θ.
- **Proposition B2** — The Falconer heritability function h²(σ) = κ²σ²/(2κ²σ² + ε²) is strictly monotone increasing on (0, ∞), proved algebraically via cross-multiplication.
- **Corollary B3** — Sufficient condition for sign reversal: when ρ < −(κ_γ σ_γ)/(2√2 κ_θ σ_θ), calibration deflates ŝ_θ < σ_θ.

To build (requires Lean 4 toolchain):

```bash
cd lean
lake build
```

## Authors

[TBD]

## License

MIT
