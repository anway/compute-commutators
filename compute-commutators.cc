#include "compute-commutators.h"

namespace compute_commutators {

typedef compute_commutators::ComputeCommutators::single_coeffs single_coeffs;

ComputeCommutators::ComputeCommutators(int n) : num_orbitals(n) {} 

void ComputeCommutators::AddInitialTerms()
{
  for (int p = 0; p <= num_orbitals; ++p) {
    for (int q = 0; q <= num_orbitals; ++q) {
      term curr_term {p, -1 * q}, curr_coeff_term {p, -1 * q};

      single_coeffs curr_coeff;
      std::set<std::vector<int> > prod_of_coeffs;
      prod_of_coeffs.insert(curr_coeff_term);
      curr_coeff.product_of_coeffs = prod_of_coeffs;
      std::vector<single_coeffs> sum_of_coeffs;
      sum_of_coeffs.push_back(curr_coeff);

      initial_terms_to_coefficients[curr_term] = sum_of_coeffs;

      for (r = 0; r <= num_orbitals; ++r) {
        for (s = 0; s <= num_orbitals; ++s) {
        }
      }
    }
  } 
}

void ComputeCommutators::InterleaveTerms()
{
}

void CalculateTrotterError()
{
}

}  // namespace compute-commutators
