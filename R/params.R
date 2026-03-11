# ===========================================================================
# params.R — Single source of truth for ALL simulation parameters
# ===========================================================================

PARAMS <- list(
  # --- Master seed ---
  MASTER_SEED = 66829L,

  # --- Simulation size ---
  N_PAIRS = 50000L,

  # --- Gompertz-Makeham ---
  B_GOMP = 0.085,
  MU_THETA = log(5e-5),
  M_EX_HIST = 0.004,
  T_MAX = 120,
  CUTOFF_AGE = 15,
  N_CALIB_ITER = 40,

  # --- Targets ---
  TARGET_H2 = 0.50,

  # --- Extrinsic heterogeneity ---
  SIGMA_GAMMA_DEF = 0.40,
  RHO_DEF = 0.4,

  # --- MGG model (Danish cohort, Shenhar Table S1) ---
  MGG_A0 = 1e-5,
  MGG_B0 = 0.115,
  MGG_C = 30,
  MGG_MU_Q = 1.0,
  MGG_CV_Q = 0.27,

  # --- SR model (Danish cohort, Shenhar Table S1) ---
  SR_ETA = 0.66,
  SR_BETA = 62.96,
  SR_EPSILON = 51.83,
  SR_KAPPA = 0.50,
  SR_MU_XC = 17.0,
  SR_CV_XC = 0.21,
  SR_DT = 1 / 52,

  # --- Shenhar's published targets (Danish, Herskind 1996) ---
  SHENHAR_R_MZ_DANISH = 0.21,
  SHENHAR_R_DZ_DANISH = 0.10,

  # --- Vanishing m_ex for control B ---
  M_EX_TINY = 0.0002,

  # --- Hamilton model ---
  HAM_E_WEIGHT = 0.61,
  HAM_EXT_SLOPE = 0.62,
  HAM_EXT_SHIFT = 0.70,
  HAM_TARGET_H2 = 0.45,

  # --- MC replication counts ---
  N_REPS_MAIN = 50L,     # main arms replication (50 seeds)
  N_REPS_SWEEP = 20L,    # sensitivity sweep replication (20 seeds/point)

  # --- Inner parallelism cap (for mclapply inside {crew} workers) ---
  # SR sweeps use mclapply for grid points; cap to avoid oversubscription
  # with {crew}'s 10 workers. 4 cores × ~2 concurrent SR targets = 8, safe on 14.
  INNER_MC_CORES = 4L
)
