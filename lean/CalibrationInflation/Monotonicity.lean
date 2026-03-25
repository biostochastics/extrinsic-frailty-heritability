/-
  CalibrationInflation.Monotonicity
  Proposition B2 (Falconer h² monotonicity) and Lemma (r_MZ calibration monotonicity).

  Both use the same functional form f(σ) = κ²σ²/(2κ²σ² + ε²).
  Proved algebraically via cross-multiplication — no calculus needed.

  Manuscript reference: Appendix B, §B.5 "Monotonicity of extrapolated h²(σ)"
-/
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum
import Mathlib.Data.Real.Basic

namespace CalibrationInflation

/-- The MZ correlation / Falconer heritability function:
    f(σ) = κ²σ² / (2κ²σ² + ε²)
    Same functional form for r_MZ^fit(sh_θ) and h²(σ). -/
noncomputable def corrFun (κ ε σ : ℝ) : ℝ :=
  (κ ^ 2 * σ ^ 2) / (2 * κ ^ 2 * σ ^ 2 + ε ^ 2)

/-- Denominator is always strictly positive. -/
theorem corrFun_denom_pos (κ ε σ : ℝ) (_hκ : κ > 0) (hε : ε > 0) :
    2 * κ ^ 2 * σ ^ 2 + ε ^ 2 > 0 := by
  have : ε ^ 2 > 0 := sq_pos_of_pos hε
  nlinarith [sq_nonneg κ, sq_nonneg σ]

/-- Denominator is nonzero (needed by field_simp). -/
theorem corrFun_denom_ne_zero (κ ε σ : ℝ) (hκ : κ > 0) (hε : ε > 0) :
    2 * κ ^ 2 * σ ^ 2 + ε ^ 2 ≠ 0 :=
  ne_of_gt (corrFun_denom_pos κ ε σ hκ hε)

/-- corrFun is strictly less than 1/2. -/
theorem corrFun_lt_half (κ ε σ : ℝ) (hκ : κ > 0) (hε : ε > 0) :
    corrFun κ ε σ < 1 / 2 := by
  unfold corrFun
  have hd := corrFun_denom_pos κ ε σ hκ hε
  rw [div_lt_div_iff₀ hd (by norm_num : (0:ℝ) < 2)]
  nlinarith [sq_pos_of_pos hε]

/-- corrFun is non-negative. -/
theorem corrFun_nonneg (κ ε σ : ℝ) (hκ : κ > 0) (hε : ε > 0) :
    corrFun κ ε σ ≥ 0 := by
  unfold corrFun
  apply div_nonneg
  · nlinarith [sq_nonneg κ, sq_nonneg σ]
  · linarith [corrFun_denom_pos κ ε σ hκ hε]

/-- corrFun is strictly positive when σ > 0. -/
theorem corrFun_pos (κ ε σ : ℝ) (hκ : κ > 0) (hε : ε > 0) (hσ : σ > 0) :
    corrFun κ ε σ > 0 := by
  unfold corrFun
  apply div_pos
  · have := sq_pos_of_pos hκ; have := sq_pos_of_pos hσ; nlinarith
  · exact corrFun_denom_pos κ ε σ hκ hε

/-- **Proposition B2 + Lemma (algebraic proof).**
    corrFun is strictly monotone increasing on (0, ∞).

    Proof: for 0 < x < y, cross-multiply denominators (both positive)
    and reduce to κ²ε²x² < κ²ε²y², which follows from x² < y². -/
theorem corrFun_strictMonoOn (κ ε : ℝ) (hκ : κ > 0) (hε : ε > 0) :
    StrictMonoOn (corrFun κ ε ·) (Set.Ioi 0) := by
  intro x hx y hy hxy
  simp only [Set.mem_Ioi] at hx hy
  unfold corrFun
  have hdx := corrFun_denom_pos κ ε x hκ hε
  have hdy := corrFun_denom_pos κ ε y hκ hε
  rw [div_lt_div_iff₀ hdx hdy]
  have hε2 : ε ^ 2 > 0 := sq_pos_of_pos hε
  have hκ2 : κ ^ 2 > 0 := sq_pos_of_pos hκ
  have hκε : κ ^ 2 * ε ^ 2 > 0 := mul_pos hκ2 hε2
  have hxy2 : x ^ 2 < y ^ 2 := by nlinarith [sq_nonneg (y - x)]
  nlinarith [mul_lt_mul_of_pos_left hxy2 hκε]

/-- Prop B2: Falconer heritability h²(σ) is strictly increasing for σ > 0. -/
theorem prop_B2_h2_strictMono (κ ε : ℝ) (hκ : κ > 0) (hε : ε > 0) :
    StrictMonoOn (corrFun κ ε ·) (Set.Ioi 0) :=
  corrFun_strictMonoOn κ ε hκ hε

/-- Lemma: r_MZ^fit is strictly increasing with range ⊆ (0, 1/2). -/
theorem lemma_rMZ_strictMono (κ ε : ℝ) (hκ : κ > 0) (hε : ε > 0) :
    StrictMonoOn (corrFun κ ε ·) (Set.Ioi 0) :=
  corrFun_strictMonoOn κ ε hκ hε

theorem lemma_rMZ_range (κ ε σ : ℝ) (hκ : κ > 0) (hε : ε > 0) (hσ : σ > 0) :
    corrFun κ ε σ ∈ Set.Ioo 0 (1 / 2) :=
  ⟨corrFun_pos κ ε σ hκ hε hσ, corrFun_lt_half κ ε σ hκ hε⟩

end CalibrationInflation
