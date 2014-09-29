#include "compute-commutators.h"

int main(int argc, char* argv[])
{
  if (argc < 2 || argc > 3) {
    printf("Usage: ./compute-commutators num_orbitals for verbose version or"
        "./compute-commutators num_orbitals F for non-verbose version.\n");
    exit(1);
  } 
  bool verbose;
  if (argc == 2) {
    verbose = true;
  } else if (*argv[2] == 'F' || *argv[2] == 'f') {
    verbose = false;
  } else {
    printf("Usage: ./compute-commutators num_orbitals for verbose version or"
        "./compute-commutators num_orbitals F for non-verbose version.\n");
    exit(2);
  }
  compute_commutators::ComputeCommutators compute_commutators =
      compute_commutators::ComputeCommutators(std::atoi(argv[1]),
      verbose);
  compute_commutators.AddInitialTerms();
  compute_commutators.InterleaveTerms();
  compute_commutators.CalculateTrotterError();
  return 0;
}
