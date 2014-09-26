#ifndef COMPUTE_COMMUTATORS_UTIL_H_
#define COMPUTE_COMMUTATORS_UTIL_H_

#include <map>
#include <set>
#include <string>
#include <vector>

namespace compute_commutators_util {

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
  static std::vector<single_coeffs> GetInitialSumCoeffs(
      std::vector<int> curr_coeff_term);
};

class TermsToCoeffsMap {
 public:
  void AddNormalForm(term curr_term, std::vector<single_coeffs> curr_coeff);
 private:
  std::map<term, std::vector<single_coeffs> > terms_to_coefficients;
};

}  // namespace compute_commutators_util

#endif
