# Qupled

Qupled is a Python package designed for calculating the properties of quantum plasmas using the dielectric formalism. By combining a straightforward Python interface with the speed of C++, it allows for efficient and accurate computations of quantum plasma properties.

<p align="center">
  <img src="examples/readme/qupled_animation_light.svg#gh-light-mode-only" width="100%">
  <img src="examples/readme/qupled_animation_dark.svg#gh-dark-mode-only" width="100%">
<p>

## Status
[![Build & Test](https://github.com/fedluc/qupled/actions/workflows/build-and-test.yml/badge.svg)](https://github.com/fedluc/qupled/actions/workflows/build-and-test.yml)
[![Code Formatting](https://github.com/fedluc/qupled/actions/workflows/formatting.yml/badge.svg)](https://github.com/fedluc/qupled/actions/workflows/formatting.yml)
![](https://readthedocs.org/projects/qupled/badge/?version=latest&style=flat)
 
## Running 

After [installation](https://qupled.readthedocs.io/en/latest/introduction.html#dependencies) qupled can be used as a regular Python package

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

<ol>
  <li>
    <a href="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.109.125134">Tolias, P., Lucco Castello, F., Kalkavouras, F., &#38; Dornheim, T. (2024). Revisiting the Vashishta-Singwi dielectric scheme for the warm dense uniform electron fluid. <i>Phys. Rev. B</i>, <i>109</i>(12)</a>
  </li>
  <li>
    <a href="https://pubs.aip.org/aip/jcp/article/158/14/141102/2877795/Quantum-version-of-the-integral-equation-theory">Tolias, P., Lucco Castello, F., &#38; Dornheim, T. (2023). Quantum version of the integral equation theory-based dielectric scheme for strongly coupled electron liquids. <i>The Journal of Chemical Physics</i>, <i>158</i>(14)</a>
  </li>
  <li>
    <a href="https://pubs.aip.org/aip/jcp/article/155/13/134115/353165/Integral-equation-theory-based-dielectric-scheme">Tolias, P., Lucco Castello, F., &#38; Dornheim, T. (2021). Integral equation theory based dielectric scheme for strongly coupled electron liquids. <i>The Journal of Chemical Physics</i>, <i>155</i>(13)</a>
  </li>
</ol>
