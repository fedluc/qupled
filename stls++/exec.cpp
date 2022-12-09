#include <omp.h>
#include <iostream>
#include "input.hpp"
#include "chemicalpotential.hpp"

int main(int argc, char** argv) {

  // Read input 
  if (argc != 2) {
    return 1;
  }
  Input in;
  try {
    in.readInput(argv[1]);
    in.print();
  }
  catch (const runtime_error& err) {
    cout << err.what() << endl;
    return 1;
  }

  // Compute chemical potential
  shared_ptr<StaticInput> statIn = in.getStaticInput();
  vector<double> guess = statIn->getChemicalPotentialGuess();
  ChemicalPotential mu(in.getDegeneracy());
  mu.compute(guess);
  cout << "Chemical potential = " << mu.get() << endl;
  return 0;
}
