#include "compute-commutators.h"
#include "compute-commutators-util.h"

#include <ctime>

namespace compute_commutators {

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
          // Only add term if it preserves spin parity
          if ((r % 2) + (s % 2) == (p % 2) + (q % 2)) {
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
  } 
  initial_terms_to_coefficients.RemoveComplexConjugates();
  // Put all the keys into initial_terms
  for (auto it = initial_terms_to_coefficients.Begin();
      it != initial_terms_to_coefficients.End(); ++it) {
    initial_terms.insert(it->first);
  }
  fprintf(stderr, "Initial terms size: %lu\n", initial_terms.size());
}

void ComputeCommutators::AddTermToInterleavedOrder(const term& curr_term) {
  if (initial_terms_to_coefficients.HasTerm(curr_term)) {
    auto it = initial_terms.find(curr_term);
    if (it != initial_terms.end()) { 
      interleaved_order.push_back(curr_term);
      initial_terms.erase(it);
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
    for (int q = 1; q <= p; ++q) {
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
    fprintf(stderr, "Order of terms in Trotter series:\n");
    for (const auto& curr_term : interleaved_order) {
      fprintf(stderr, "Term: ");
      compute_commutators_util::ComputeCommutatorsUtil::PrintIndices(stderr,
          curr_term);
      fprintf(stderr, "\nCoeff: ");
      compute_commutators_util::ComputeCommutatorsUtil::PrintSumOfCoeffs(
          stderr, initial_terms_to_coefficients.At(curr_term));
      fprintf(stderr, "\n");
    }
  }
  fprintf(stderr, "Number of terms in interleaved order: %lu\n",
      interleaved_order.size());
}

std::pair<std::vector<term>, all_coeff> ComputeCommutators
    ::GetTermForTrotter(const int& index) {
  term curr_term = interleaved_order[index];
  term curr_conjugate = compute_commutators_util::ComputeCommutatorsUtil
      ::GetConjugate(curr_term);
  all_coeff curr_coeff = initial_terms_to_coefficients.At(curr_term);
  std::vector<term> term_and_conjugate;
  if (!curr_conjugate.empty()) {
    // return both term and conjugate if conjugate and term differ
    term_and_conjugate.push_back(curr_term);
    term_and_conjugate.push_back(curr_conjugate);
  } else {
    // otherwise return just term
    term_and_conjugate.push_back(curr_term);
  }
  return std::make_pair(term_and_conjugate, curr_coeff);
}

void ComputeCommutators::CalculateTrotterError() {
  int num_commutators = 0, num_terms = interleaved_order.size();
  for (int i = 0; i < num_terms; ++i) {
    num_commutators += i * (i + 1);  
  }
  int one_percent = num_commutators / 100.0;
  fprintf(stderr, "There are %d possible commutators.\n", num_commutators);
  int counter = 0;
  clock_t start = clock();

  // Loop over the possible combinations a <= b, c < b of
  // (1 / 12) * [A * (1 - delta(A, B)/2), [B, C]] = ...
  // ((1 - delta(A, B) / 2) / 12) * (ABC - ACB - BCA + CBA)

  // Loop over B.
  for (int b = 0; b < num_terms; ++b) {
    const auto& B_term_coeff = GetTermForTrotter(b);
    std::vector<term> B_terms = B_term_coeff.first;
    all_coeff B_coeff = B_term_coeff.second;

    // Loop over A.
    for (int a = 0; a <= b; ++a) {
      const auto& A_term_coeff = GetTermForTrotter(a);
      std::vector<term> A_terms = A_term_coeff.first;
      all_coeff A_coeff = A_term_coeff.second;

      // Loop over C.
      for (int c = 0; c < b; ++c) {
        const auto& C_term_coeff = GetTermForTrotter(c);
        std::vector<term> C_terms = C_term_coeff.first;
        all_coeff C_coeff = C_term_coeff.second;
       
        // Increment counter.
        counter += 1;

        // Calculate the coefficient (multiplying the 3 terms). 
        all_coeff multiplied_coeffs = compute_commutators_util
            ::ComputeCommutatorsUtil::MultiplySumOfCoeffs(A_coeff, B_coeff);
        multiplied_coeffs = compute_commutators_util::ComputeCommutatorsUtil
            ::MultiplySumOfCoeffs(multiplied_coeffs, C_coeff);
        // To handle the scale (each coefficient needs to be multipled by 1 / 24
        // if a == b, or 1 / 12 if not): we multiply the coefficient by 2 if
        // a is not equal to b, and divide everything by 24 in the very end.
        // This allows us to keep coefficients as ints. 
        if (a == b) {
          for (auto it = multiplied_coeffs.begin();
              it != multiplied_coeffs.end(); ++it) {
            it->second *= 2;
          }
        }

        // Compute commutators.
        for (const term& A : A_terms) {
          for (const term& B : B_terms) {
            for (const term& C : C_terms) {
              if (!compute_commutators_util::ComputeCommutatorsUtil
                  ::TriviallyCommutes(A, B, C)) {
                // ABC term
                final_terms_to_coefficients.AddNormalForm(
                    compute_commutators_util::ComputeCommutatorsUtil::
                    ConcatenateThreeTerms(A, B, C), multiplied_coeffs);
                // CBA term
                final_terms_to_coefficients.AddNormalForm(
                    compute_commutators_util::ComputeCommutatorsUtil::
                    ConcatenateThreeTerms(C, B, A), multiplied_coeffs);
                // Flip sign on coefficients for ACB and BCA terms.
                for (auto it = multiplied_coeffs.begin();
                    it != multiplied_coeffs.end(); ++it) {
                  it->second *= -1;
                }
                // ACB term
                final_terms_to_coefficients.AddNormalForm(
                    compute_commutators_util::ComputeCommutatorsUtil::
                    ConcatenateThreeTerms(A, C, B), multiplied_coeffs);
                // BCA term
                final_terms_to_coefficients.AddNormalForm(
                    compute_commutators_util::ComputeCommutatorsUtil::
                    ConcatenateThreeTerms(B, C, A), multiplied_coeffs);
              }
            }
          }
        }
      
        // Report progress.
        if (one_percent != 0 && counter % one_percent == 0) {
          int percent_complete = counter / one_percent;
          int elapsed = double(clock() - start) / CLOCKS_PER_SEC; 
          int rate = elapsed / percent_complete; 
          int eta = rate * (100 - percent_complete); 
          // Get current date and time.
          time_t rawtime;
          time(&rawtime);
          struct tm *timeinfo = localtime(&rawtime);
          fprintf(stderr, "%sComputation %d%% complete. Approximately "
              "%d minute(s) remaining.\n", asctime(timeinfo),
              percent_complete, eta / 60);
        }
      } 
    }
  } 
}

void ComputeCommutators::PrintFinalResults(FILE* output) {
  for (auto it = final_terms_to_coefficients.Begin();
      it != final_terms_to_coefficients.End(); ++it) {
    compute_commutators_util::ComputeCommutatorsUtil::PrintIndices(output,
        it->first);
    fprintf(output, "\n");
    compute_commutators_util::ComputeCommutatorsUtil::PrintSumOfCoeffs(output,
        it->second);
    fprintf(output, "\n\n");
  }
}

}  // namespace compute-commutators
