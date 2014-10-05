#ifndef COMPUTE_COMMUTATORS_H_
#define COMPUTE_COMMUTATORS_H_

#include "compute-commutators-util.h"

#include <set>
#include <string>
#include <vector>

namespace compute_commutators {

typedef std::vector<int> term;
typedef compute_commutators_util::all_coeff all_coeff;
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
  void AddTermToInterleavedOrder(const term& curr_term);
  // Calculate the Trotter error using the terms arranged in the order from
  // InterleaveTerms, and with the complex conjugates added back in.
  void CalculateTrotterError();
  // Helper for CalculateTrotterError to return term, its conjugate, its coeff.
  std::pair<std::vector<term>, all_coeff > GetTermForTrotter(
      const int& index);
  void PrintFinalResults(FILE* output);
 private:
  int num_orbitals;
  TermsToCoeffsMap initial_terms_to_coefficients;
  // initial_terms is an intermediate set we will use to construct
  // interleaved_order
  std::set<term> initial_terms;
  std::vector<term> interleaved_order;
  TermsToCoeffsMap final_terms_to_coefficients;
  bool verbose;
};

} //  namespace compute_commutators

#endif
