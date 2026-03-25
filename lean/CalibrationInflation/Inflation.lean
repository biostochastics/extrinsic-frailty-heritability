/-
  CalibrationInflation.Inflation
  Proposition B1: inflation under non-negative pleiotropy.

  The central theorem. Under feasibility + non-negative pleiotropy (ρ ≥ 0),
  the calibrated intrinsic dispersion exceeds the true value: sh_θ > σ_θ.

  Manuscript reference: Appendix B, Proposition B1 (line 1214)
-/
import CalibrationInflation.Params
import CalibrationInflation.Moments
import CalibrationInflation.Calibration
import Mathlib.Tactic.FieldSimp

namespace CalibrationInflation

/-- The core algebraic inequality underlying Proposition B1.
    Under non-negative pleiotropy, the numerator of (V̂ - V_θ) is positive:
    (V_γ + C_ρ) · σ_ε² + V_θ · C_ρ > 0 -/
theorem inflation_numerator_pos (p : Params) (hρ : p.ρ ≥ 0) :
    (p.V_γ + p.C_ρ) * p.σ_ε ^ 2 + p.V_θ * p.C_ρ > 0 := by
  have hVγ := p.V_γ_pos
  have hCρ := p.C_ρ_nonneg hρ
  have hσε := p.σ_ε_sq_pos
  have hVθ := p.V_θ_pos
  nlinarith

/-- Inflation of calibrated variance:
    κ_θ² sh_θ² > κ_θ² σ_θ²
    under non-negative pleiotropy and feasibility. -/
theorem calibrated_variance_inflated (p : Params)
    (hρ : p.ρ ≥ 0)
    (h_feas : p.Feasible)
    (sh_θ : ℝ)
    (h_calib : p.CalibrationEq sh_θ) :
    calibratedVariance p.κ_θ sh_θ > p.V_θ := by
  have h_exp := calibrated_variance_expanded p h_feas sh_θ h_calib
  rw [h_exp]
  have h_denom_pos : p.σ_ε ^ 2 - p.C_ρ > 0 := by
    rw [← p.V_true_minus_2C]; exact h_feas
  rw [gt_iff_lt, lt_div_iff₀ h_denom_pos]
  nlinarith [inflation_numerator_pos p hρ]

/-- **Proposition B1 (local inflation under first-order approximation).**

    Omission of heritable extrinsic frailty (σ_γ > 0) with non-negative
    pleiotropy (ρ ≥ 0) inflates the calibrated intrinsic dispersion: sh_θ > σ_θ. -/
theorem prop_B1_inflation (p : Params)
    (hρ : p.ρ ≥ 0)
    (h_feas : p.Feasible)
    (sh_θ : ℝ)
    (h_calib : p.CalibrationEq sh_θ)
    (h_sh_pos : sh_θ > 0) :
    sh_θ > p.σ_θ := by
  have h_var := calibrated_variance_inflated p hρ h_feas sh_θ h_calib
  unfold calibratedVariance Params.V_θ at h_var
  have hκ2 : p.κ_θ ^ 2 > 0 := sq_pos_of_pos p.hκ_θ
  have h_sq : sh_θ ^ 2 > p.σ_θ ^ 2 := by nlinarith
  nlinarith [sq_nonneg (sh_θ - p.σ_θ)]

end CalibrationInflation
