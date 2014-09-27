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

term ComputeCommutatorsUtil::GetConjugate(term curr_term) {
  term conjugate;
  for (int index : curr_term) {
    conjugate.push_back(-1 * index);
  }
  if (conjugate == curr_term) {
    conjugate.clear();
  } 
  return conjugate;
}

void TermsToCoeffsMap::AddNormalForm(term curr_term,
    std::vector<single_coeffs> curr_coeff) {
  for (std::vector<int>::iterator it = curr_term.begin(); it != curr_term.end();
      ++it) {
    for (std::vector<int>::iterator rit = it; rit != curr_term.begin();
        --rit) {
      const int right = *rit, left = *(rit - 1);
      if (right > left) {
        *(rit - 1) = right;
        *rit = left;
        if (left == -1 * right) {
          term exchange_term;
          exchange_term.insert(exchange_term.end(), curr_term.begin(), rit - 1);
          exchange_term.insert(exchange_term.end(), rit + 1, curr_term.end());
          AddNormalForm(exchange_term, curr_coeff);
        }
        for (single_coeffs one_coeff_term : curr_coeff) {
          one_coeff_term.integer_multiplier *= -1;
        }
      } else {
        break;
      }
    }
  }
  std::set<int> no_repeats(curr_term.begin(), curr_term.end());
  if (no_repeats.size() == curr_term.size()) {
    if (terms_to_coefficients.find(curr_term) == terms_to_coefficients.end()) {
      terms_to_coefficients[curr_term] = curr_coeff;
    } else {
      terms_to_coefficients[curr_term].insert(
          terms_to_coefficients[curr_term].end(), curr_coeff.begin(),
          curr_coeff.end());
    }
  }
}

void TermsToCoeffsMap::RemoveComplexConjugates() {
  for (const auto& term_to_coeff : terms_to_coefficients) {
    term conjugate = compute_commutators_util::ComputeCommutatorsUtil
        ::GetConjugate(term_to_coeff.first);
    if (!conjugate.empty()) {
      auto it = terms_to_coefficients.find(conjugate);
      if (it != terms_to_coefficients.end()) {
        terms_to_coefficients.erase(it);
      }
    }    
  }
}

}  // namespace compute_commutators_util
