#ifndef LAUNCHER_HPP
#define LAUNCHER_HPP

class RpaInput;

class MyOutput {
public:
  MyOutput(int rank_) : rank(rank_) { ; }
  int getRank() const { return rank; }
private:
  int rank;
};

Rpa mpiExampleFunction(const RpaInput& in);

#endif
