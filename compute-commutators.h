#ifndef COMPUTE_COMMUTATORS_H_
#define COMPUTE_COMMUTATORS_H_

#include "compute-commutators-util.h"

#include <string>
#include <vector>

namespace compute_commutators {

typedef std::vector<int> term;
typedef compute_commutators_util::single_coeffs single_coeffs;
typedef compute_commutators_util::TermsToCoeffsMap TermsToCoeffsMap;

class ComputeCommutators {
 public:
  ComputeCommutators(int n, bool verbose);
  // Add initial terms in Hamiltonian (pq and pqrs) along with symbolic
  // representation of coefficients to the terms to coefficients maps with 
  // terms normal ordered (so swapping takes place as necessary), and then
  // remove complex conjugates from map.
  void AddInitialTerms();
  // Prepare the order of terms for the final calculation of Trotter error.
  void InterleaveTerms();
  // Calculate the Trotter error using the terms arranged in the order from
  // InterleaveTerms, and with the complex conjugates added back in.
  void CalculateTrotterError();
 private:
  int num_orbitals;
  TermsToCoeffsMap initial_terms_to_coefficients;
  std::vector<term> interleaved_order;
  TermsToCoeffsMap final_terms_to_coefficients;
  bool verbose;
};

} //  namespace compute_commutators

#endif
