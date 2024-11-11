# Qupled

Qupled is a Python package designed for calculating the properties of quantum plasmas using the dielectric formalism. By combining a straightforward Python interface with the speed of C++, it allows for efficient and accurate computations of quantum plasma properties.

<p align="center">
  <img src="examples/readme/qupled_animation_light.svg#gh-light-mode-only" width="100%">
  <img src="examples/readme/qupled_animation_dark.svg#gh-dark-mode-only" width="100%">
<p>

## Status

[![Build & Test (Linux)](https://github.com/fedluc/qupled/actions/workflows/build_and_test_ubuntu.yml/badge.svg)](https://github.com/fedluc/qupled/actions/workflows/build_and_test_ubuntu.yml)
[![Build & Test (Linux-MPI)](https://github.com/fedluc/qupled/actions/workflows/build_and_test_ubuntu_mpi.yml/badge.svg)](https://github.com/fedluc/qupled/actions/workflows/build_and_test_ubuntu_mpi.yml)
[![Build & Test (macOS)](https://github.com/fedluc/qupled/actions/workflows/build_and_test_macos.yml/badge.svg)](https://github.com/fedluc/qupled/actions/workflows/build_and_test_macos.yml)
[![Build & Test (macOS-MPI)](https://github.com/fedluc/qupled/actions/workflows/build_and_test_macos_mpi.yml/badge.svg)](https://github.com/fedluc/qupled/actions/workflows/build_and_test_macos_mpi.yml)
[![Code Formatting](https://github.com/fedluc/qupled/actions/workflows/formatting.yml/badge.svg)](https://github.com/fedluc/qupled/actions/workflows/formatting.yml)

## Dependencies

Before trying to build, test or run the code one should make sure that [these dependencies](https://qupled.readthedocs.io/en/latest/introduction.html#dependencies) are satisfied.
 
## Building & running

Qupled can be compiled with `cmake`, tested with `pytest` and installed with the following procedure

```bash
git clone https://github.com/fedluc/qupled.git
cd qupled
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
pytest tests
cmake --install .
```
After installation Qupled can be used as a regular Python package

```python
# Solve the stls dielectric scheme for coupling = 10 and degeneracy 1.0
from qupled.classic import Stls
inputs = Stls.Input(10.0, 1.0)
Stls().compute(inputs)
```

## Documentation

More detailed information on the package together with a list of examples is available in the [documentation](http://qupled.readthedocs.io/)

## Publications

Qupled has been used in the following publications:

``` bibtex
@article{tolias2021integral,
  title={Integral equation theory based dielectric scheme for strongly coupled electron liquids},
  author={Tolias, Panagiotis and Lucco Castello, F and Dornheim, Tobias},
  journal={The Journal of Chemical Physics},
  volume={155},
  number={13},
  year={2021},
  publisher={AIP Publishing}
}

@article{tolias2023quantum,
  title={Quantum version of the integral equation theory-based dielectric scheme for strongly coupled electron liquids},
  author={Tolias, Panagiotis and Lucco Castello, Federico and Dornheim, Tobias},
  journal={The Journal of Chemical Physics},
  volume={158},
  number={14},
  year={2023},
  publisher={AIP Publishing}
}

@article{PhysRevB.109.125134,
  title = {Revisiting the Vashishta-Singwi dielectric scheme for the warm dense uniform electron fluid},
  author = {Tolias, Panagiotis and Lucco Castello, Federico and Kalkavouras, Fotios and Dornheim, Tobias},
  journal = {Phys. Rev. B},
  volume = {109},
  issue = {12},
  pages = {125134},
  numpages = {22},
  year = {2024},
  month = {Mar},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.109.125134},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.109.125134}
}


```
