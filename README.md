# Qupled: quantum plasmas and dielectric formalism

Platform | CI Status
---------|:---------
Linux      | [![OSX Build & Test Status](https://github.com/fedluc/qupled/actions/workflows/build_and_test_ubuntu.yml/badge.svg)](https://github.com/fedluc/qupled/actions/workflows/build_and_test_ubuntu.yml)
OSX      | [![OSX Build & Test Status](https://github.com/fedluc/qupled/actions/workflows/build_and_test_macos.yml/badge.svg)](https://github.com/fedluc/qupled/actions/workflows/build_and_test_macos.yml)


Qupled is a python package that can be used to compute the properties of quantum plasmas via the dielectric formalism. The plasma properties can be computed to arbitrary precision by leveraging on a simple Python interface combined with the speed of C++

<p align="center">
 <img src="examples/readme/qupled_animation_light.gif#gh-light-mode-only" width="800" height="533">
 <img src="examples/readme/qupled_animation_dark.gif#gh-dark-mode-only" width="800" height="533">
<p>
 
## Building & running

Qupled can be compiled with `cmake`, tested with `pytest` and installed with the following procedure

```bash
git clone git@github.com:fedluc/qupled.git
cd qupled
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
python3 -m pip install pandas tables pytest pytest-mock 
pytest tests
cmake --install .
```
After installation Qupled can be used as a regular Python package

```python
import qupled.classic as qpc
import qupled.quantum as qpq
# Solve the stls dielectric scheme for coupling = 10 and degeneracy 1.0
qpc.Stls(10.0, 1.0).compute()
# Solve the qstls dielectric scheme for coupling = 5 and degeneracy 2.0
qpq.Qstls(5.0, 2.0).compute()
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
