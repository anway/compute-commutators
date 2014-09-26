#include "compute-commutators.h"
#include "compute-commutators-util.h"

namespace compute_commutators {

typedef compute_commutators_util::single_coeffs single_coeffs;

ComputeCommutators::ComputeCommutators(int n) : num_orbitals(n) {} 

void ComputeCommutators::AddInitialTerms() {
  for (int p = 0; p <= num_orbitals; ++p) {
    for (int q = 0; q <= num_orbitals; ++q) {
      term curr_term;
      curr_term.push_back(p);
      curr_term.push_back(-1 * q);
      term curr_coeff_term(curr_term);

      initial_terms_to_coefficients.AddNormalForm(curr_term,
          compute_commutators_util::ComputeCommutatorsUtil::GetInitialSumCoeffs(
          curr_coeff_term));

      for (int r = 0; r <= num_orbitals; ++r) {
        for (int s = 0; s <= num_orbitals; ++s) {
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
}

void ComputeCommutators::InterleaveTerms() {
}

void CalculateTrotterError() {
}

}  // namespace compute-commutators
