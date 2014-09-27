#include "compute-commutators.h"
#include "compute-commutators-util.h"

namespace compute_commutators {

typedef compute_commutators_util::single_coeffs single_coeffs;

ComputeCommutators::ComputeCommutators(int n, bool verbose) : num_orbitals(n),
    verbose(verbose) {} 

void ComputeCommutators::AddInitialTerms() {
  for (int p = 1; p <= num_orbitals; ++p) {
    for (int q = 1; q <= num_orbitals; ++q) {
      // pq term.
      term curr_term;
      curr_term.push_back(p);
      curr_term.push_back(-1 * q);
      term curr_coeff_term(curr_term);

      initial_terms_to_coefficients.AddNormalForm(curr_term,
          compute_commutators_util::ComputeCommutatorsUtil::GetInitialSumCoeffs(
          curr_coeff_term));

      for (int r = 1; r <= num_orbitals; ++r) {
        for (int s = 1; s <= num_orbitals; ++s) {
          // pqrs term.
          term curr_term;
          curr_term.push_back(p);
          curr_term.push_back(q);
          curr_term.push_back(-1 * r);
          curr_term.push_back(-1 * s);
          term curr_coeff_term(curr_term);

          initial_terms_to_coefficients.AddNormalForm(curr_term,
              compute_commutators_util::ComputeCommutatorsUtil
              ::GetInitialSumCoeffs(curr_coeff_term));
        }
      }
    }
  } 
  initial_terms_to_coefficients.RemoveComplexConjugates();
}

void ComputeCommutators::InterleaveTerms() {
  // First add Hpp terms.
  for (int p = 1; p <= num_orbitals; ++p) {
    term curr_term;
    curr_term.push_back(p);
    curr_term.push_back(-1 * p);

    if (initial_terms_to_coefficients.HasTerm(curr_term)) {
      interleaved_order.push_back(curr_term);
    }
  } 
  // Now add Hpqqp terms.
}

void CalculateTrotterError() {
}

}  // namespace compute-commutators
