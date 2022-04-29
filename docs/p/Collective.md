---
title: Collective Effects
---

A collective effects subpackage `pyat/at/collective` allows to model impedance driven collective effects and perform multi-particle tracking. Presently only single bunch effects are included. It is possible to build pyAT with MPI to perform computing intensive simulations on a cluster.

### Wake object

The `Wake` provides an interface to create an object containing the wake field information. It is then passed to the lattice element that is used for tracking. A `Wake` is defined by its `s` coordinate and wake componenents: transverse dipole, transverse quadrupole and longitudinal.
