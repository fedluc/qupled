#include <iostream>
#include <mpi.h>
#include "input.hpp"
#include "rpa.hpp"
#include "launcher.hpp"



Rpa mpiExampleFunction(const RpaInput& in) {
  int rank, size;
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // MyOutput output(-1); // Initialize with a sentinel value
  // if (rank == 0) {
  //   // Only rank 0 performs the operation
  //   // Add your MPI-related logic here
  //   std::cout << "Hello from MPI process " << rank << "coupling is  " << in.getCoupling() << std::endl;
  //   // Update the MyOutput instance with the rank
  //   output = MyOutput(rank);
  // }
  // Broadcast the result from rank 0 to all ranks
  //MPI_Bcast(&output, sizeof(MyOutput), MPI_BYTE, 0, MPI_COMM_WORLD);
  Rpa output(in);
  MPI_Finalize();
  return output;
}
