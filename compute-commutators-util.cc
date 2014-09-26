#include "compute-commutators.h"

namespace compute_commutators_util {

typedef compute_commutators_util::single_coeffs single_coeffs;

std::vector<single_coeffs> ComputeCommutatorsUtil::GetInitialSumCoeffs(
    std::vector<int> curr_coeff_term) {
  single_coeffs curr_coeff;
  std::set<std::vector<int> > prod_of_coeffs;
  prod_of_coeffs.insert(curr_coeff_term);
  curr_coeff.product_of_coeffs = prod_of_coeffs;
  std::vector<single_coeffs> sum_of_coeffs;
  sum_of_coeffs.push_back(curr_coeff);

  return sum_of_coeffs;
}

void TermsToCoeffsMap::AddNormalForm(term curr_term,
    std::vector<single_coeffs> curr_coeff) {
}

}  // namespace compute_commutators_util
