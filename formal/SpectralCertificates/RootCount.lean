import Mathlib.Algebra.Polynomial.Roots

/-!
# Kernel-checked scalar branch obstruction

This file formalizes the polynomial root-count step used in the exact
invariant-line obstruction.  Once tangential full-spark geometry supplies
more than `N` distinct pass nodes on which a scalar branch equals one, a
degree-at-most-`N` branch is identically one.
-/

namespace SpectralCertificates

open Polynomial

theorem scalar_branch_eq_one_of_many_pass_nodes
    (p : ℝ[X]) (N : ℕ) (nodes : Finset ℝ)
    (degree_le : p.natDegree ≤ N)
    (many_nodes : N < nodes.card)
    (calibrated : ∀ x ∈ nodes, p.eval x = 1) :
    p = 1 := by
  apply Polynomial.eq_of_natDegree_lt_card_of_eval_eq' p 1 nodes calibrated
  simp only [Polynomial.natDegree_one, Nat.max_zero]
  omega

end SpectralCertificates
