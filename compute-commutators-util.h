#ifndef COMPUTE_COMMUTATORS_UTIL_H_
#define COMPUTE_COMMUTATORS_UTIL_H_

#include <map>
#include <set>
#include <string>
#include <vector>

namespace compute_commutators_util {

// Each term is represented by a list of indices, where a positive index
// represents a raising operator and a negative index represents a lowering
// operator.
typedef std::vector<int> term;

// single_coeff represents a term of the form integer_number * h_1 * h_2 * ...
// where the h_i are symbolic representations of h_{pq}, h_{pqrs}
// We will then represent a sum of single_coeff with a vector of single_coeffs
struct single_coeffs {
  int integer_multiplier = 1;
  std::set<std::vector<int> > product_of_coeffs;
};

class ComputeCommutatorsUtil {
 public:
  // A helper that returns the coefficient corresponding to the initial term.
  // The coefficient should have the same form as the initial term (i.e., same
  // list of indices [p,-q] or [p,q,-r,-s]) since we have not done any swaps
  // or multiplications yet at this point.
  static std::vector<single_coeffs> GetInitialSumCoeffs(
      std::vector<int> curr_coeff_term);
  // Returns the conjugate of a term in normal order, or an empty list is a term
  // is its own conjugate.
  static term GetConjugate(term curr_term);
};

class TermsToCoeffsMap {
 public:
  // Adds a term and coefficient to the map, with the term in normal order
  // after necessary swapping.
  void AddNormalForm(term curr_term, std::vector<single_coeffs> curr_coeff);
  // Removes complex conjugates from the map.
  void RemoveComplexConjugates();
  bool HasTerm(term curr_term);
 private:
  std::map<term, std::vector<single_coeffs> > terms_to_coefficients;
};

}  // namespace compute_commutators_util

#endif
