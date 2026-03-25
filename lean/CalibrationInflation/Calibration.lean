/-
  CalibrationInflation.Calibration
  Calibration equation solution existence and expanded form.

  Manuscript reference: Appendix B, §B.4 "Calibration equation and inflation condition"
-/
import CalibrationInflation.Params
import CalibrationInflation.Moments
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Tactic.FieldSimp

namespace CalibrationInflation

/-- The RHS of the calibration equation is positive under feasibility with C_true > 0. -/
theorem calibration_rhs_pos (p : Params) (h_feas : p.Feasible)
    (h_C_pos : p.C_true > 0) :
    (p.C_true * p.σ_ε ^ 2) / (p.V_true - 2 * p.C_true) > 0 := by
  apply div_pos
  · have := p.σ_ε_sq_pos; nlinarith
  · exact h_feas

/-- The calibration equation has a positive solution.
    Explicitly: sh_θ = √(C_true · σ_ε² / (κ_θ² · (V_true - 2·C_true))). -/
theorem calibration_exists (p : Params) (h_feas : p.Feasible)
    (h_C_pos : p.C_true > 0) :
    ∃ sh_θ : ℝ, sh_θ > 0 ∧ p.CalibrationEq sh_θ := by
  set rhs := (p.C_true * p.σ_ε ^ 2) / (p.κ_θ ^ 2 * (p.V_true - 2 * p.C_true)) with rhs_def
  have h_rhs_pos : rhs > 0 := by
    rw [rhs_def]
    apply div_pos
    · have := p.σ_ε_sq_pos; nlinarith
    · exact mul_pos (sq_pos_of_pos p.hκ_θ) h_feas
  exact ⟨Real.sqrt rhs, Real.sqrt_pos.mpr h_rhs_pos, by
    unfold Params.CalibrationEq calibratedVariance
    rw [Real.sq_sqrt (le_of_lt h_rhs_pos), rhs_def]
    have hκ2 : p.κ_θ ^ 2 ≠ 0 := ne_of_gt (sq_pos_of_pos p.hκ_θ)
    have hd : p.V_true - 2 * p.C_true ≠ 0 := p.feasible_ne_zero h_feas
    have hκ_ne : p.κ_θ ≠ 0 := p.κ_θ_ne_zero
    field_simp [hκ2, hd, hκ_ne]⟩

/-- Calibrated variance in expanded form (manuscript eq. @calib-expanded):
    κ_θ² sh_θ² = (V_θ + V_γ + C_ρ) · σ_ε² / (σ_ε² - C_ρ) -/
theorem calibrated_variance_expanded (p : Params) (_h_feas : p.Feasible) (sh_θ : ℝ)
    (h_calib : p.CalibrationEq sh_θ) :
    calibratedVariance p.κ_θ sh_θ =
      (p.V_θ + p.V_γ + p.C_ρ) * p.σ_ε ^ 2 / (p.σ_ε ^ 2 - p.C_ρ) := by
  unfold Params.CalibrationEq at h_calib
  rw [h_calib, p.V_true_minus_2C]
  unfold Params.C_true
  rfl

end CalibrationInflation
