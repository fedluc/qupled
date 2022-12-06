#include <omp.h>
#include <input.hpp>
#include <iostream>

int main(int argc, char** argv) {
  if (argc != 2) {
    return 1;
  }
  Input input;
  try {
    input.readInput(argv[1]);
    input.print();
  }
  catch (const runtime_error& err) {
    cout << err.what() << endl;
    return 1;
  }
  return 0;
}
