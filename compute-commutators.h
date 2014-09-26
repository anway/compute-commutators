#ifndef COMPUTE_COMMUTATORS_H_
#define COMPUTE_COMMUTATORS_H_

#include "compute-commutators-util.h"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace compute_commutators {

typedef std::vector<int> term;
typedef compute_commutators_util::single_coeffs single_coeffs;

class ComputeCommutators {
 public:
  explicit ComputeCommutators(int n);
  void AddInitialTerms();
  void InterleaveTerms();
  void CalculateTrotterError();
  /*
  // single_coeff represents a term of the form integer_number * h_1 * h_2 *...
  // where the h_i are symbolic representations of h_{pq}, h_{pqrs}
  // We will then represent a sum of single_coeff with a vector of single_coeffs
  struct single_coeffs {
    int integer_multiplier = 1;
    std::set<std::vector<int> > product_of_coeffs;
  };
  */
 private:
  int num_orbitals;
  std::map<term, std::vector<single_coeffs> >
      initial_terms_to_coefficients;
  std::vector<term> interleaved_order;
  std::map<term, std::vector<single_coeffs> >
      final_terms_to_coefficients;
};

} //  namespace compute_commutators

#endif
