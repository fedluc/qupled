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
      //in.print();
    }
    catch (const runtime_error& err) {
      cerr << err.what() << endl;
      return 1;
    }
  }


  // // Compute slts scheme
  // Stls stls(in);
  // stls.compute();

  // Compute qstls scheme
  Qstls qstls(in);
  qstls.compute();
  return 0;
  
}
