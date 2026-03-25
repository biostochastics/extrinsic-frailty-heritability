/-
  CalibrationInflation.Params
  Parameter structure for the calibration inflation model.

  Manuscript reference: Appendix B, §B.1 "Setup: true DGP and fitted model"
  Source: Kornilov, "Omitted familial extrinsic risk can inflate
          inferred intrinsic lifespan heritability"
-/
import Mathlib.Analysis.SpecialFunctions.Sqrt
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum

namespace CalibrationInflation

/-- Parameters for the calibration inflation model.
    All scale parameters strictly positive. Pleiotropy correlation ρ ∈ [-1, 1]. -/
structure Params where
  /-- Intrinsic sensitivity coefficient (κ_θ in manuscript) -/
  κ_θ : ℝ
  /-- Extrinsic sensitivity coefficient (κ_γ in manuscript) -/
  κ_γ : ℝ
  /-- True intrinsic frailty scale -/
  σ_θ : ℝ
  /-- Extrinsic frailty scale -/
  σ_γ : ℝ
  /-- Residual variance scale (equal-noise approximation: σ_{ε,T} ≈ σ_{ε,F}) -/
  σ_ε : ℝ
  /-- Pleiotropy correlation: Corr(G, γ^gen) -/
  ρ : ℝ
  hκ_θ : κ_θ > 0
  hκ_γ : κ_γ > 0
  hσ_θ : σ_θ > 0
  hσ_γ : σ_γ > 0
  hσ_ε : σ_ε > 0
  hρ_le : ρ ≤ 1
  hρ_ge : ρ ≥ -1

-- Derived quantities

/-- Intrinsic variance component: V^T_θ = κ_θ² σ_θ² -/
def Params.V_θ (p : Params) : ℝ := p.κ_θ ^ 2 * p.σ_θ ^ 2

/-- Extrinsic variance component: V^T_γ = (1/2) κ_γ² σ_γ² -/
noncomputable def Params.V_γ (p : Params) : ℝ := (1 / 2) * p.κ_γ ^ 2 * p.σ_γ ^ 2

/-- Pleiotropy cross-term: C_ρ = √2 ρ κ_θ κ_γ σ_θ σ_γ -/
noncomputable def Params.C_ρ (p : Params) : ℝ := Real.sqrt 2 * p.ρ * p.κ_θ * p.κ_γ * p.σ_θ * p.σ_γ

-- √2 helper lemmas

theorem sqrt2_pos : (0 : ℝ) < Real.sqrt 2 := by positivity

theorem sqrt2_sq : (Real.sqrt 2) ^ 2 = 2 :=
  Real.sq_sqrt (by norm_num : (2 : ℝ) ≥ 0)

-- Positivity lemmas (structure fields pulled into context via `have`)

theorem Params.V_θ_pos (p : Params) : p.V_θ > 0 := by
  have := p.hκ_θ; have := p.hσ_θ
  unfold V_θ; positivity

theorem Params.V_γ_pos (p : Params) : p.V_γ > 0 := by
  have := p.hκ_γ; have := p.hσ_γ
  unfold V_γ; positivity

theorem Params.C_ρ_nonneg (p : Params) (hρ : p.ρ ≥ 0) : p.C_ρ ≥ 0 := by
  have := p.hκ_θ; have := p.hκ_γ; have := p.hσ_θ; have := p.hσ_γ
  unfold C_ρ; positivity

theorem Params.σ_ε_sq_pos (p : Params) : p.σ_ε ^ 2 > 0 := by
  have := p.hσ_ε; positivity

-- Nonzero lemmas (needed by field_simp)

theorem Params.κ_θ_ne_zero (p : Params) : p.κ_θ ≠ 0 := ne_of_gt p.hκ_θ
theorem Params.σ_θ_ne_zero (p : Params) : p.σ_θ ≠ 0 := ne_of_gt p.hσ_θ
theorem Params.σ_ε_ne_zero (p : Params) : p.σ_ε ≠ 0 := ne_of_gt p.hσ_ε

end CalibrationInflation
