Introduction
============

qupled is a package that can be used to compute the properties of quantum one component
plasmas via theoretical approaches based on the dielectric formalism. The theoretical
approaches which can be solved with qupled include:

  * The classical STLS scheme as discussed by `Tanaka and Ichimaru <https://journals.jps.jp/doi/abs/10.1143/JPSJ.55.2278>`_
  * The classical STLS-HNC scheme as discussed by `Tanaka <https://pubs.aip.org/aip/jcp/article/145/21/214104/196066/Correlational-and-thermodynamic-properties-of>`_
  * The classical STLS-IET scheme as discussed by `Tolias and collaborators <https://pubs.aip.org/aip/jcp/article/155/13/134115/353165/Integral-equation-theory-based-dielectric-scheme>`_
  * The quantum STLS (QSTLS) scheme as discussed by `Schweng and Böhm <https://journals.aps.org/prb/abstract/10.1103/PhysRevB.48.2037>`_ 
  * The quantum QSTLS-IET scheme as discussed by `Tolias and others <https://pubs.aip.org/aip/jcp/article/158/14/141102/2877795/Quantum-version-of-the-integral-equation-theory>`_

Limitations
===========

Ground state (zero temperature) calculations are available only for the classical schemes (STLS, STLS-HNC and STLS-IET).

Units
=====

All the calculations are performed in normalized units. The wave vectors are normalized to the Fermi wave-vector and the frequencies are normalized to 2πEF/h. Here EF is the Fermi energy and h is Planck's constant.
