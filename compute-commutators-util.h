#ifndef COMPUTE_COMMUTATORS_UTIL_H_
#define COMPUTE_COMMUTATORS_UTIL_H_

namespace compute_commutators_util {

// single_coeff represents a term of the form integer_number * h_1 * h_2 * ...
// where the h_i are symbolic representations of h_{pq}, h_{pqrs}
// We will then represent a sum of single_coeff with a vector of single_coeffs
struct single_coeffs {
  int integer_multiplier = 1;
  std::set<std::vector<int> > product_of_coeffs;
};

class ComputeCommutatorsUtil {
 public:
  static vector<

}  // namespace compute_commutators_util
