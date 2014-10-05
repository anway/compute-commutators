#include "compute-commutators-util.h"

namespace compute_commutators_util {

typedef compute_commutators_util::single_coeff single_coeff;

all_coeff ComputeCommutatorsUtil::GetInitialSumCoeffs(term curr_coeff_term) {
  // Return a sum of coeffs that has just one coeff in the sum (the
  // one corresponding to the initial term).
  single_coeff one_coeff;
  all_coeff sum_coeffs_map;
  // If a pq term, check for pq = qp symmetry. If a pqrs term, check for
  // pqrs = qprs = pqsr = qpsr symmetry 
  if (curr_coeff_term.size() == 2) {
    // swap so that we consistently have [p, q] with p < q
    if (curr_coeff_term[1] < curr_coeff_term[0]) {
      int temp = curr_coeff_term[1];
      curr_coeff_term[1] = curr_coeff_term[0];
      curr_coeff_term[0] = temp;
    }
  } else if (curr_coeff_term.size() == 4) {
    // swap so that we consistently have [p, q, r, s] with p < q, r < s
    if (curr_coeff_term[1] < curr_coeff_term[0]) {
      int temp = curr_coeff_term[1];
      curr_coeff_term[1] = curr_coeff_term[0];
      curr_coeff_term[0] = temp;
    }
    if (curr_coeff_term[3] < curr_coeff_term[2]) {
      int temp = curr_coeff_term[3];
      curr_coeff_term[3] = curr_coeff_term[2];
      curr_coeff_term[2] = temp;
    }
  }
  // Just one term in product.
  one_coeff.insert(curr_coeff_term);
  // Insert one coefficient with multiplier 1 into map.
  sum_coeffs_map[one_coeff] = 1;

  return sum_coeffs_map;
}

term ComputeCommutatorsUtil::GetConjugate(const term& curr_term) {
  term conjugate;
  for (auto rit = curr_term.rbegin(); rit != curr_term.rend(); ++rit) {
    conjugate.push_back(-1 * (*rit));
  }
  if (conjugate == curr_term) {
    // Return nothing if a term is equal to its own conjugate.
    conjugate.clear();
  } 
  return conjugate;
}

void ComputeCommutatorsUtil::PrintIndices(FILE* output, const term& curr_term) {
  fprintf(output, "[ ");
  for (int index : curr_term) {
    fprintf(output, "%d ", index);
  }
  fprintf(output, "] ");
}

void ComputeCommutatorsUtil::PrintSumOfCoeffs(FILE* output,
    const all_coeff& sum_of_coeffs) {
  for (auto it = sum_of_coeffs.begin(); it != sum_of_coeffs.end();
      ++it) {
    fprintf(output, "%d * ", it->second);
    for (const auto& one_term : it->first) {
      PrintIndices(output, one_term);
    } 
    if (it != --sum_of_coeffs.end()) {
      fprintf(output, " + ");
    }
  }
}

all_coeff ComputeCommutatorsUtil::MultiplySumOfCoeffs(
    const all_coeff& first, const all_coeff& second) {
  // result will carry result of multiplication
  all_coeff result;
  for (auto first_term_it = first.begin(); first_term_it != first.end();
      ++first_term_it) {
    for (auto second_term_it = second.begin(); second_term_it != second.end();
        ++second_term_it) {
      const auto& first_term = first_term_it->first;
      const auto& second_term = second_term_it->first;
      // first multiply just the symbolic terms (without coefficient)
      std::multiset<std::vector<int> > new_term;
      new_term.insert(first_term.begin(), first_term.end());
      new_term.insert(second_term.begin(), second_term.end());    
      // check if resulting symbolic term has already been seen in this
      // multiplication
      auto it_result = result.find(new_term);
      // if not, push new symbolic term and coefficient
      if (it_result == result.end()) {
        result[new_term] = first_term_it->second * second_term_it->second;
      } else {  // if so, just add to coefficient of existing term
        it_result->second += first_term_it->second * second_term_it->second;  
      }
    }
  }
  return result;
}

bool ComputeCommutatorsUtil::TriviallyCommutes(const term& first_term,
    const term& second_term, const term& third_term) {
  // First check if B and C have indices in common
  std::set<int> set_one(second_term.begin(), second_term.end());
  std::set<int> set_two(third_term.begin(), third_term.end());
  std::vector<int> intersection;
  std::set_intersection(set_one.begin(), set_one.end(), set_two.begin(),
      set_two.end(), std::back_inserter(intersection));
  if (!intersection.empty()) {
    return false;
  }
  intersection.clear();
  // Now check if A and BC have indices in common
  std::copy(third_term.begin(), third_term.end(), std::inserter(
      set_one, set_one.end()));
  std::set<int> set_three(first_term.begin(), first_term.end());
  std::set_intersection(set_one.begin(), set_one.end(), set_three.begin(),
      set_three.end(), std::back_inserter(intersection));
  if (!intersection.empty()) {
    return false;
  } else {
    return true;
  }
}

term ComputeCommutatorsUtil::ConcatenateThreeTerms(const term& first_term,
    const term& second_term, const term& third_term) {
  term result;
  result.insert(result.end(), first_term.begin(), first_term.end());
  result.insert(result.end(), second_term.begin(), second_term.end());
  result.insert(result.end(), third_term.begin(), third_term.end());
  return result;
}

void TermsToCoeffsMap::AddNormalForm(term curr_term, all_coeff curr_coeff) {
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
          if (!exchange_term.empty()) {
            AddNormalForm(exchange_term, curr_coeff);
          }
        }
        // Flip the sign of the coefficient because we swapped.
        for (auto it_coeff = curr_coeff.begin(); it_coeff != curr_coeff.end();
           ++it_coeff) {
          it_coeff->second *= -1;
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
    } else {  // Term already in map; just add (sum) existing coefficients.
      for (auto it_coeff = curr_coeff.begin(); it_coeff != curr_coeff.end();
         ++it_coeff) {
        // If the multiplier is 0, we can ignore it
        if (it_coeff->second != 0) {
          // Check if we're just combining existing terms or adding new term.
          auto full_map_it = terms_to_coefficients[curr_term].find(
              it_coeff->first);
          // If not, push new term and coefficient.
          if (full_map_it == terms_to_coefficients[curr_term].end()) {
            (terms_to_coefficients[curr_term])[it_coeff->first] =
                it_coeff->second; 
          } else {
            // Add integer multipliers if term already in map.
            full_map_it->second += it_coeff->second;
            if (full_map_it->second == 0) {
              terms_to_coefficients[curr_term].erase(full_map_it);
            }
          }
        }
      }
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

bool TermsToCoeffsMap::HasTerm(const term& curr_term) {
  if (terms_to_coefficients.find(curr_term) != terms_to_coefficients.end()) {
    return true;
  } else {
    return false;
  }
}

std::map<term, all_coeff>::iterator TermsToCoeffsMap
    ::Begin() {
  return terms_to_coefficients.begin();
}

std::map<term, all_coeff>::iterator TermsToCoeffsMap
    ::End() {
  return terms_to_coefficients.end();
}

all_coeff TermsToCoeffsMap::At(const term& key) {
  return terms_to_coefficients.at(key);
}

}  // namespace compute_commutators_util
