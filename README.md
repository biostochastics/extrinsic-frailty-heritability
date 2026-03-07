# Omitted extrinsic frailty inflates intrinsic heritability in calibrate-then-extrapolate mortality models

Simulation pipeline demonstrating that **model misspecification** (omitted heritable extrinsic frailty) biases Shenhar et al.'s intrinsic heritability estimate upward.

## Summary

[Shenhar et al. (2026)](https://doi.org/10.1126/science.adz1187) estimate that the heritability of human lifespan rises to ~50% after removing extrinsic mortality, using parametric mortality models calibrated to Scandinavian twin data. [Hamilton (2026)](https://doi.org/10.64898/2026.02.26.26347172) argues that infection mortality is heritable and demonstrates collider bias from concordant-survivor conditioning. However, Shenhar's procedure does not remove individuals — it parametrically sets extrinsic hazard to zero.

We identify the bias mechanism specific to Shenhar's calibrate-then-extrapolate procedure: if individuals vary in heritable susceptibility to extrinsic mortality and this heterogeneity is not modeled, the single-parameter calibration absorbs extrinsic genetic variance into the intrinsic parameter — omitted-variable bias at the estimation step.

### Key results (Gompertz–Makeham)

| Condition | σ\_θ,fit | h² | Bias (pp) |
|---|---|---|---|
| Oracle (true intrinsic) | 1.259 | 0.510 | — |
| Arm 1: Correctly specified | 1.220 | 0.472 | −3.7 |
| **Arm 2: Misspecification** (σ\_γ=0.40, ρ=0.4) | **1.537** | **0.586** | **+7.6** |
| Control A: Non-heritable extrinsic | 1.259 | 0.478 | −3.1 |
| Control B: Vanishing m\_ex | 1.267 | 0.500 | −0.9 |
| Pleiotropy isolation (ρ=0) | 1.376 | 0.528 | +1.9 |
| Irrelevant trait (ζ) | 1.264 | 0.502 | −0.7 |
| Oracle fix (correct model) | 1.219 | 0.476 | −3.3 |
| Cutoff = 0 (no age restriction) | 1.335 | 0.514 | +3.6 |
| Arm 3: Hamilton conditioning | N/A | 0.427 | −8.2 |

The bias decomposes into ~3.7 pp of downward attenuation (estimand mismatch) and +11.4 pp of upward inflation (variance absorption), netting +7.6 pp.

### Model-generality validation

| Model | Oracle h² | Arm 1 h² | Arm 1 bias | Arm 2 h² | Arm 2 bias |
|---|---|---|---|---|---|
| Gompertz–Makeham | 0.510 | 0.472 | −3.7 pp | 0.586 | **+7.6 pp** |
| MGG (Shenhar) | 0.464 | 0.459 | −0.5 pp | 0.556 | **+9.2 pp** |
| SR (Shenhar) | 0.577 | 0.569 | −0.8 pp | 0.655 | **+7.8 pp** |

### Anchored sensitivity regime

Parameters anchored to independent immunogenetic evidence: σ\_γ from twin heritability of infection mortality ([Obel et al. 2010](https://doi.org/10.1093/aje/kwq037); h²≈0.40) and ρ from LDSC genetic correlations between infection severity and longevity ([Qiu et al. 2025](https://doi.org/10.1186/s12967-024-05932-y); r\_g≈−0.35).

| Region | Bias range (pp) | Mean bias (pp) |
|---|---|---|
| Externally informed regime (σ\_γ∈[0.30,0.65], ρ∈[0.20,0.50]) | +2.8 to +18.9 | +9 |
| Conservative point estimate (σ\_γ=0.65, ρ=0.35) | — | +15.7 |

## Repository structure

```
├── _targets.R              # Pipeline definition (183 targets, static branching)
├── R/
│   ├── params.R            # Single source of truth for all parameters
│   ├── sim_core.R          # Lambert W GM solver, hybrid MGG solver, SR Euler-Maruyama
│   ├── sim_twins.R         # Twin parameter generation, sim_twin_h2()
│   ├── calibrate.R         # Stochastic bisection (CRN), oracle, joint, multi-target
│   ├── utils.R             # Sweep cells, Hamilton arm, model analysis, control sweeps
│   ├── anchors.R           # σ_γ bridge, tetrachoric, m_ex split, negative ρ
│   ├── diagnostics.R       # Variance decomposition, joint r_MZ/r_DZ, α/β estimation
│   ├── plots.R             # All figure-generating functions (16+ plots)
│   └── export.R            # CSV/JSON export
├── src/
│   └── sim_lifespan_sr.cpp # Rcpp/C++ SR solver (Xoroshiro128+ RNG, ~4.6× speedup)
├── tables/                 # Pipeline outputs (CSV/JSON)
│   ├── scalars.json        # All scalar results
│   ├── summary.csv         # Main summary table
│   ├── models.csv          # Cross-model validation
│   ├── sweep.csv           # Full 81-cell sensitivity grid
│   ├── anchored.csv        # 56-cell anchored grid
│   └── controls.csv        # Control experiments
└── figures/                # All generated figures (PNG)
```

## Reproducing

### Clone and install

```bash
git clone https://github.com/biostochastics/extrinsic-frailty-heritability.git
cd extrinsic-frailty-heritability
```

**R** (≥ 4.4):
```r
install.packages(c("targets", "tarchetypes", "crew",
                    "MASS", "ggplot2", "patchwork", "lamW",
                    "ggpubr", "ggthemes", "viridis", "mvtnorm",
                    "scales", "Rcpp"))
```

### Run the pipeline

```bash
# Full pipeline (~20-30 min on Apple M4 Pro, 10 parallel workers via crew)
Rscript -e 'targets::tar_make()'

# Check pipeline status
Rscript -e 'targets::tar_visnetwork()'

# Inspect a specific target
Rscript -e 'targets::tar_read(arm2_corr)'
```

The pipeline uses `{crew}` for parallel execution with 10 local workers. It defines 183 targets organized in a DAG with `tar_map` static branching for the 81-cell sensitivity sweep and 56-cell anchored sweep. A master seed (66829) propagates deterministic sub-seeds to every target for full reproducibility.

### Pipeline architecture

```
Stage 0:  Calibration         → sigma_theta_true, oracle
Stage 1:  Main arms           → Arm 1, Arm 2, Controls A/B, Arm 3, Oracle fix
Stage 2:  Sensitivity sweep   → 81 cells (tar_map static branching)
Stage 3:  Anchored sweep      → 56 cells (tar_map static branching)
Stage 4:  Variance decomp     → α, β coefficients
Stage 5:  Model validation    → MGG + SR arms
Stage 6:  Additional controls → Dose-response, pleiotropy isolation, irrelevant trait
Stage 7:  Diagnostics         → Joint r_MZ/r_DZ, σ_γ bridge, m_ex split, negative ρ
Stage 8:  Calibration UQ      → Bootstrap uncertainty
Stage 9:  Tables              → CSV/JSON export
Stage 10: Figures             → All plots
```

## Key technical notes

- **GM solver**: Closed-form via Lambert W function (`lamW` package), with bisection fallback for edge cases.
- **MGG parameterization**: We use `a = a₀^q` (Strehler-Mildvan compensation), not `a = a₀^(1/q)` as stated in Shenhar's supplement. The latter does not produce hazard curve crossing and yields inconsistent twin correlations.
- **SR solver**: Euler-Maruyama at weekly resolution (Δt = 1/52), compiled to C++ via Rcpp with Xoroshiro128+ RNG for ~4.6× speedup over vectorized R.
- **Calibration**: 40-iteration stochastic bisection with common random numbers (CRN) for variance reduction.

## Author

Sergey A. Kornilov, Biostochastics, Seattle, WA

## Acknowledgements

GPT-5.2-Pro (OpenAI) and Gemini-3-Pro (Google) were used during development to review the simulation design, verify that the calibration procedure parallels Shenhar et al.'s methodology, and provide structural feedback. MiniMax-M2.5 (MiniMax), GLM-5 (Zhipu AI), GPT-5.2-Pro, and Gemini-3-Pro were used to review and assist in deriving the formal proofs. Claude Opus 4.6 (Anthropic) was used for code generation, drafting assistance, and proof typesetting. All code, analytical results, and final text were reviewed and revised by the author.

## License

MIT
