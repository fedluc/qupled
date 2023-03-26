#include <omp.h>
#include <iostream>
#include "input.hpp"
#include "stls.hpp"
#include "qstls.hpp"

int main(int argc, char** argv) {
  // Read input
  Input in;
  if (argc < 2) {
    cout << "Running with all default values" << endl;
  }
  else {
    try {
      in.readInput(argv[1]);
    }
    catch (const runtime_error& err) {
      cerr << err.what() << endl;
      return 1;
    }
  }
  // Set number of omp threads
  omp_set_num_threads(in.getNThreads());
  // Compute scheme
  try {
    if (in.isClassic()) {
      Stls stls(in);
      stls.compute();
    }
    else {
      Qstls qstls(in);
      qstls.compute();
    }
  }
  catch (const runtime_error& err) {
    cerr << err.what() << endl;
    return 1;
  }
  // Return success
  return 0;
}
