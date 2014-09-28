#include "compute-commutators.h"

namespace compute_commutators_util {

typedef compute_commutators_util::single_coeffs single_coeffs;

std::vector<single_coeffs> ComputeCommutatorsUtil::GetInitialSumCoeffs(
    std::vector<int> curr_coeff_term) {
  // Return a sum of single_coeffs that has just one coeff in the sum (the
  // one corresponding to the initial term).
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
    // Return nothing if a term is equal to its own conjugate.
    conjugate.clear();
  } 
  return conjugate;
}

void ComputeCommutatorsUtil::PrintIndices(term curr_term) {
  printf("[ ");
  for (int index : curr_term) {
    printf("%d ", index);
  }
  printf("] ");
}

void ComputeCommutatorsUtil::PrintSumOfCoeffs(std::vector<single_coeffs>
    sum_of_coeffs) {
  for (auto it = sum_of_coeffs.begin(); it != sum_of_coeffs.end();
      ++it) {
    printf("%d * ", (*it).integer_multiplier);
    for (const auto& prod_of_terms : (*it).product_of_coeffs) {
      PrintIndices(prod_of_terms);
    }  
    if (it != sum_of_coeffs.end() - 1) {
      printf(" + ");
    }
  }
}

void TermsToCoeffsMap::AddNormalForm(term curr_term,
    std::vector<single_coeffs> curr_coeff) {
  for (std::vector<int>::iterator it = curr_term.begin(); it != curr_term.end();
      ++it) {
    for (std::vector<int>::iterator rit = it; rit != curr_term.begin();
        --rit) {
      const int right = *rit, left = *(rit - 1);
      if (right > left) {
        // Swap the two indices.
        *(rit - 1) = right;
        *rit = left;
        if (left == -1 * right) {
          // If additionally we are swapping two indices that act on the same
          // tensor factor, add an additional term without the two indices.
          // (fermionic commutation relations)
          term exchange_term;
          exchange_term.insert(exchange_term.end(), curr_term.begin(), rit - 1);
          exchange_term.insert(exchange_term.end(), rit + 1, curr_term.end());
          AddNormalForm(exchange_term, curr_coeff);
        }
        // Flip the sign of the coefficient because we swapped.
        for (single_coeffs one_coeff_term : curr_coeff) {
          one_coeff_term.integer_multiplier *= -1;
        }
      } else {  // No more terms out of order, so stop swapping.
        break;
      }
    }
  }
  std::set<int> no_repeats(curr_term.begin(), curr_term.end());
  // If an index is repeated, we have two raising or two lowering operators, so
  // don't add to map.
  if (no_repeats.size() == curr_term.size()) {
    if (terms_to_coefficients.find(curr_term) == terms_to_coefficients.end()) {
      // Term not in map already.
      terms_to_coefficients[curr_term] = curr_coeff;
    } else {  // Term already in map; just add coefficients.
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
    if (!conjugate.empty()) {  // If the term is not its own conjugate
      auto it = terms_to_coefficients.find(conjugate);
      // Remove its conjugate from map.
      if (it != terms_to_coefficients.end()) {
        terms_to_coefficients.erase(it);
      }
    }    
  }
}

bool TermsToCoeffsMap::HasTerm(term curr_term) {
  if (terms_to_coefficients.find(curr_term) != terms_to_coefficients.end()) {
    return true;
  } else {
    return false;
  }
}

std::map<term, std::vector<single_coeffs> >::iterator TermsToCoeffsMap
    ::Begin() {
  return terms_to_coefficients.begin();
}

std::map<term, std::vector<single_coeffs> >::iterator TermsToCoeffsMap
    ::End() {
  return terms_to_coefficients.end();
}

std::vector<single_coeffs> TermsToCoeffsMap::At(term key) {
  return terms_to_coefficients.at(key);
}

}  // namespace compute_commutators_util
