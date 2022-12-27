#include <omp.h>
#include <iostream>
#include "input.hpp"
#include "stls.hpp"

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
    cerr << err.what() << endl;
    return 1;
  }

  // Compute slts scheme
  Stls stls(in);
  stls.compute();
  return 0;
}
