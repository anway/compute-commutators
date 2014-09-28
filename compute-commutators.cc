#include "compute-commutators.h"
#include "compute-commutators-util.h"

namespace compute_commutators {

typedef compute_commutators_util::single_coeffs single_coeffs;

ComputeCommutators::ComputeCommutators(int n, bool verbose) : num_orbitals(n),
    verbose(verbose) {} 

void ComputeCommutators::AddInitialTerms() {
  for (int p = 1; p <= num_orbitals; ++p) {
    for (int q = 1; q <= num_orbitals; ++q) {
      // pq term.
      term curr_term;
      curr_term.push_back(p);
      curr_term.push_back(-1 * q);
      term curr_coeff_term(curr_term);

      initial_terms_to_coefficients.AddNormalForm(curr_term,
          compute_commutators_util::ComputeCommutatorsUtil::GetInitialSumCoeffs(
          curr_coeff_term));

      for (int r = 1; r <= num_orbitals; ++r) {
        for (int s = 1; s <= num_orbitals; ++s) {
          // pqrs term.
          term curr_term;
          curr_term.push_back(p);
          curr_term.push_back(q);
          curr_term.push_back(-1 * r);
          curr_term.push_back(-1 * s);
          term curr_coeff_term(curr_term);

          initial_terms_to_coefficients.AddNormalForm(curr_term,
              compute_commutators_util::ComputeCommutatorsUtil
              ::GetInitialSumCoeffs(curr_coeff_term));
        }
      }
    }
  } 
  initial_terms_to_coefficients.RemoveComplexConjugates();
  // Put all the keys into initial_terms
  for (auto it = initial_terms_to_coefficients.Begin();
      it != initial_terms_to_coefficients.End(); ++it) {
    initial_terms.insert(it->first);
  }
}

void ComputeCommutators::AddTermToInterleavedOrder(term curr_term) {
  if (initial_terms_to_coefficients.HasTerm(curr_term)) {
    auto it = initial_terms.find(curr_term);
    if (it != initial_terms.end()) { 
      interleaved_order.push_back(curr_term);
      initial_terms.erase(initial_terms.find(curr_term));
    }
  }      
}

void ComputeCommutators::InterleaveTerms() {
  // First add Hpp terms.
  for (int p = 1; p <= num_orbitals; ++p) {
    term curr_term;
    curr_term.push_back(p);
    curr_term.push_back(-1 * p);

    AddTermToInterleavedOrder(curr_term);
  } 
  // Now add Hpqqp terms.
  for (int p = 1; p <= num_orbitals; ++p) {
    for (int q = 1; q <= p; ++q) {
      term curr_term;
      curr_term.push_back(p);
      curr_term.push_back(q);
      curr_term.push_back(-1 * q);
      curr_term.push_back(-1 * p);

      AddTermToInterleavedOrder(curr_term);
    }
  }
  // Now interleave (p, -q) and (p, r, -r, -q) terms.
  for (int p = 1; p <= num_orbitals; ++p) {
    for (int q = 1; q <= num_orbitals; ++q) {
      term curr_term;
      curr_term.push_back(p);
      curr_term.push_back(-1 * q);
  
      AddTermToInterleavedOrder(curr_term);

      for (int r = 1; r <= p && r <= q; ++r) {
        term curr_term;
        curr_term.push_back(p);
        curr_term.push_back(r);
        curr_term.push_back(-1 * r);
        curr_term.push_back(-1 * q); 

        AddTermToInterleavedOrder(curr_term);
      }
    }
  }
  // Now add pqrs terms.
  for (int p = 1; p <= num_orbitals; ++p) {
    for (int q = 1; q <= num_orbitals; ++q) {
      for (int r = 1; r <= num_orbitals; ++r) {
        for (int s = r; s <= num_orbitals; ++s) {
          term curr_term;
          curr_term.push_back(p);
          curr_term.push_back(q);
          curr_term.push_back(-1 * r);
          curr_term.push_back(-1 * s);
           
          AddTermToInterleavedOrder(curr_term);
        }
      }
    }
  }

  if (verbose) {
    printf("Order of terms in Trotter series:\n");
    for (const auto& curr_term : initial_terms) {
      compute_commutators_util::ComputeCommutatorsUtil::PrintIndices(curr_term);
      printf(", ");
      compute_commutators_util::ComputeCommutatorsUtil::PrintSumOfCoeffs(
          initial_terms_to_coefficients.At(curr_term));
      printf("\n");
    }
  }
}

std::pair<term, std::vector<single_coeffs> > ComputeCommutators
    ::GetTermForTrotter(int index) {
  term curr_term = interleaved_order[index];
  term curr_conjugate = compute_commutators_util::ComputeCommutatorsUtil
      ::GetConjugate(curr_term);
  std::vector<single_coeffs> curr_coeff =
      initial_terms_to_coefficients.At(curr_term);
  if (!curr_conjugate.empty()) {
    curr_term.insert(curr_term.end(), curr_conjugate.begin(),
        curr_conjugate.end());    
  }
  return std::make_pair(curr_term, curr_coeff);
}

void ComputeCommutators::CalculateTrotterError() {
  int num_commutators = 0, num_terms = interleaved_order.size();
  for (int i = 1; i < num_terms; ++i) {
    num_commutators += i * (i + 1);  
  }
  int one_percent = num_commutators / 100.0;
  printf("There are %d possible commutators.\n", num_commutators);

  // Loop over the possible combinations a <= b, c < b of
  // (1 / 12) * [A * (1 - delta(A, B)/2), [B, C]] = ...
  // ((1 - delta(A, B) / 2) / 12) * (ABC - ACB - BCA + CBA)

  // Loop over B.
  for (int b = 0; b < num_terms; ++b) {
    const auto& B_term_coeff = GetTermForTrotter(b);
    term B = B_term_coeff.first;
    std::vector<single_coeffs> B_coeff = B_term_coeff.second;

    // Loop over A.
    for (int a = 0; a <= b; ++a) {
      const auto& A_term_coeff = GetTermForTrotter(a);
      term A = A_term_coeff.first;
      std::vector<single_coeffs> A_coeff = A_term_coeff.second;

      // Loop over C.
      for (int c = 0; c < b; ++c) {
        const auto& C_term_coeff = GetTermForTrotter(c);
        term C = C_term_coeff.first;
        std::vector<single_coeffs> C_coeff = C_term_coeff.second;

        // Calculate the coefficient (multiplying the 3 terms). 
      } 
    }
  } 
}

}  // namespace compute-commutators
