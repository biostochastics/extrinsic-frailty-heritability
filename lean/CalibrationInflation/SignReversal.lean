/-
  CalibrationInflation.SignReversal
  Corollary B3: sign reversal under sufficiently negative pleiotropy.

  When ρ < -(κ_γ σ_γ)/(2√2 κ_θ σ_θ), the bias reverses: sh_θ < σ_θ.
  This is a sufficient (not necessary) condition.

  Manuscript reference: Appendix B, Corollary B3 (line 1240)
-/
import CalibrationInflation.Params
import CalibrationInflation.Moments
import CalibrationInflation.Calibration
import Mathlib.Tactic.FieldSimp

namespace CalibrationInflation

/-- Sign reversal threshold (sufficient condition):
    ρ_crit = -(κ_γ σ_γ) / (2√2 κ_θ σ_θ) -/
noncomputable def signReversalThreshold (p : Params) : ℝ :=
  -(p.κ_γ * p.σ_γ) / (2 * Real.sqrt 2 * p.κ_θ * p.σ_θ)

/-- When ρ < threshold, V_γ + C_ρ < 0. -/
theorem V_γ_plus_C_ρ_neg (p : Params) (hρ : p.ρ < signReversalThreshold p) :
    p.V_γ + p.C_ρ < 0 := by
  unfold Params.V_γ Params.C_ρ signReversalThreshold at *
  have h_denom_pos : 2 * Real.sqrt 2 * p.κ_θ * p.σ_θ > 0 := by
    have := p.hκ_θ; have := p.hσ_θ; have := sqrt2_pos; positivity
  have hρ_cleared : p.ρ * (2 * Real.sqrt 2 * p.κ_θ * p.σ_θ) <
      -(p.κ_γ * p.σ_γ) := by
    rwa [lt_div_iff₀ h_denom_pos] at hρ
  have hsq := sqrt2_sq
  have hκ_γ := p.hκ_γ; have hσ_γ := p.hσ_γ
  -- Multiply hρ_cleared by (κ_γ * σ_γ > 0) to get the key bound
  have h_prod_pos : p.κ_γ * p.σ_γ > 0 := mul_pos hκ_γ hσ_γ
  have h_mult := mul_lt_mul_of_pos_right hρ_cleared h_prod_pos
  -- h_mult: ρ*(2√2*κ_θ*σ_θ)*(κ_γ*σ_γ) < -(κ_γ*σ_γ)*(κ_γ*σ_γ)
  nlinarith [sq_nonneg (Real.sqrt 2)]

/-- **Corollary B3 (sign reversal under negative pleiotropy).**
    If ρ < -(κ_γ σ_γ)/(2√2 κ_θ σ_θ), calibration deflates: sh_θ < σ_θ.
    Sufficient (not necessary) condition for sign reversal. -/
theorem corollary_B3_deflation (p : Params)
    (hρ : p.ρ < signReversalThreshold p)
    (h_feas : p.Feasible)
    (sh_θ : ℝ)
    (h_calib : p.CalibrationEq sh_θ)
    (h_sh_pos : sh_θ > 0) :
    sh_θ < p.σ_θ := by
  -- Expanded form: V̂ = (V_θ + V_γ + C_ρ) · σ_ε² / (σ_ε² - C_ρ)
  have h_exp := calibrated_variance_expanded p h_feas sh_θ h_calib
  -- Need V̂ < V_θ, i.e. calibrated variance deflated
  suffices calibratedVariance p.κ_θ sh_θ < p.V_θ by
    unfold calibratedVariance Params.V_θ at this
    have hκ2 : p.κ_θ ^ 2 > 0 := sq_pos_of_pos p.hκ_θ
    have h_sq : sh_θ ^ 2 < p.σ_θ ^ 2 := by nlinarith
    -- sh_θ² < σ_θ² with both > 0 ⟹ sh_θ < σ_θ
    nlinarith [sq_nonneg (p.σ_θ - sh_θ), mul_pos h_sh_pos p.hσ_θ]
  rw [h_exp]
  have h_denom_pos : p.σ_ε ^ 2 - p.C_ρ > 0 := by
    rw [← p.V_true_minus_2C]; exact h_feas
  rw [div_lt_iff₀ h_denom_pos]
  -- Goal: (V_θ + V_γ + C_ρ) · σ_ε² < V_θ · (σ_ε² - C_ρ)
  -- i.e., (V_γ + C_ρ) · σ_ε² + V_θ · C_ρ < 0
  have h_neg := V_γ_plus_C_ρ_neg p hρ
  have hσε := p.σ_ε_sq_pos
  have hVθ := p.V_θ_pos
  -- C_ρ < 0 when V_γ + C_ρ < 0 (since V_γ > 0)
  have hVγ := p.V_γ_pos
  have hCρ_neg : p.C_ρ < 0 := by linarith
  -- (V_γ + C_ρ) · σ_ε² < 0 and V_θ · C_ρ < 0, so sum < 0
  nlinarith

end CalibrationInflation
