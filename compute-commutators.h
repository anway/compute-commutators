#ifndef COMPUTE_COMMUTATORS_H_
#define COMPUTE_COMMUTATORS_H_

#include "compute-commutators-util.h"

#include <set>
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
  // Helper for InterleaveTerms to add terms to interleaved_order
  void AddTermToInterleavedOrder(term curr_term);
  // Calculate the Trotter error using the terms arranged in the order from
  // InterleaveTerms, and with the complex conjugates added back in.
  void CalculateTrotterError();
  // Helper for CalculateTrotterError to prepare A/B/C term, conjugate, coeff.
  std::pair<term, std::vector<single_coeffs> > GetTermForTrotter(int index);
 private:
  int num_orbitals;
  TermsToCoeffsMap initial_terms_to_coefficients;
  std::set<term> initial_terms;
  std::vector<term> interleaved_order;
  TermsToCoeffsMap final_terms_to_coefficients;
  bool verbose;
};

} //  namespace compute_commutators

#endif
