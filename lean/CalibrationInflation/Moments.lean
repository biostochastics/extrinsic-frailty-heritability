/-
  CalibrationInflation.Moments
  Moment definitions from delta-method linearization + calibration equation.

  These are the results of the first-order linearization of log-lifespan.
  We accept them as *definitions* rather than proving from measure theory.

  Manuscript reference: Appendix B, §B.3–B.4
-/
import CalibrationInflation.Params

namespace CalibrationInflation

/-- True MZ covariance (manuscript eq. @ctrue):
    C_true = κ_θ²σ_θ² + (1/2)κ_γ²σ_γ² + √2·ρ·κ_θ·κ_γ·σ_θ·σ_γ -/
noncomputable def Params.C_true (p : Params) : ℝ :=
  p.V_θ + p.V_γ + p.C_ρ

/-- True total variance (manuscript eq. @vtrue):
    V_true = 2κ_θ²σ_θ² + κ_γ²σ_γ² + √2·ρ·κ_θ·κ_γ·σ_θ·σ_γ + σ_ε² -/
noncomputable def Params.V_true (p : Params) : ℝ :=
  2 * p.V_θ + 2 * p.V_γ + p.C_ρ + p.σ_ε ^ 2

/-- Key identity: V_true - 2·C_true = σ_ε² - C_ρ -/
theorem Params.V_true_minus_2C (p : Params) :
    p.V_true - 2 * p.C_true = p.σ_ε ^ 2 - p.C_ρ := by
  unfold V_true C_true V_θ V_γ C_ρ
  ring

/-- Feasibility condition: V_true - 2·C_true > 0.
    Required for calibration equation to have a solution. -/
noncomputable def Params.Feasible (p : Params) : Prop :=
  p.V_true - 2 * p.C_true > 0

/-- Feasibility ↔ σ_ε² > C_ρ. -/
theorem Params.feasible_iff (p : Params) :
    p.Feasible ↔ p.σ_ε ^ 2 > p.C_ρ := by
  unfold Feasible
  rw [p.V_true_minus_2C]
  constructor <;> intro h <;> linarith

/-- Feasibility holds when residual variance dominates the pleiotropy cross-term.
    Works for any ρ ∈ [-1, 1] (uses p.hρ_le internally). -/
theorem Params.feasible_of_large_residual (p : Params)
    (h_large_resid : p.σ_ε ^ 2 > Real.sqrt 2 * p.κ_θ * p.κ_γ * p.σ_θ * p.σ_γ) :
    p.Feasible := by
  rw [p.feasible_iff]
  unfold C_ρ
  have hρ_le := p.hρ_le
  have hκ_θ := p.hκ_θ; have hκ_γ := p.hκ_γ
  have hσ_θ := p.hσ_θ; have hσ_γ := p.hσ_γ
  -- ρ ≤ 1 and all factors positive → ρ * product ≤ 1 * product
  have h_prod_pos : Real.sqrt 2 * p.κ_θ * p.κ_γ * p.σ_θ * p.σ_γ > 0 := by positivity
  have h_bound : Real.sqrt 2 * p.ρ * p.κ_θ * p.κ_γ * p.σ_θ * p.σ_γ ≤
      Real.sqrt 2 * 1 * p.κ_θ * p.κ_γ * p.σ_θ * p.σ_γ := by nlinarith
  linarith

/-- Feasibility → nonzero denominator (for field_simp). -/
theorem Params.feasible_ne_zero (p : Params) (h_feas : p.Feasible) :
    p.V_true - 2 * p.C_true ≠ 0 :=
  ne_of_gt h_feas

/-- σ_ε² - C_ρ ≠ 0 when feasible (for field_simp). -/
theorem Params.σ_ε_sq_minus_C_ρ_ne_zero (p : Params) (h_feas : p.Feasible) :
    p.σ_ε ^ 2 - p.C_ρ ≠ 0 := by
  rw [← p.V_true_minus_2C]; exact p.feasible_ne_zero h_feas

/-- Calibrated variance: V̂ = κ_θ² sh_θ² -/
def calibratedVariance (κ_θ sh_θ : ℝ) : ℝ := κ_θ ^ 2 * sh_θ ^ 2

/-- Calibration equation (manuscript eq. @calib-eq):
    κ_θ² sh_θ² = C_true · σ_ε² / (V_true - 2·C_true) -/
noncomputable def Params.CalibrationEq (p : Params) (sh_θ : ℝ) : Prop :=
  calibratedVariance p.κ_θ sh_θ = (p.C_true * p.σ_ε ^ 2) / (p.V_true - 2 * p.C_true)

end CalibrationInflation
