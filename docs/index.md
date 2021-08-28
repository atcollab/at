---
title: Welcome to AT
layout: page
permalink: /
toc: false
hide_sidebar: true
---
## Introduction

Accelerator Toolbox (AT) is a collection of tools to model storage rings and beam transport lines.

With AT, it is possible to:

- **create and manipulate accelerator lattice elements**,
- **track particles through the lattice**, selecting the appropriate integrator to represent the physics
- **compute accelerator parameters and beam properties**, generating new scripts or taking advantage of the existing ones

AT is based on a 6-D modular tracking engine written in C/C++ for efficiency.
Lattice manipulation and computation of accelerator physics parameters are provided
by two interfaces:
- a [Matlab interface][1], available as a Matlab toolbox,
- a [python interface][2], available as a python package.

## Coordinate system
The 6-d phase space coordinates used in AT are as follows

$$
\begin{equation}
\vec Z = \pmatrix{x\cr \frac{p_x}{P_0} \cr y \cr \frac{p_y}{P_0} \cr \delta \cr c\tau}
\end{equation}
$$

The momenta are defined as

$$
\begin{equation}
p_x = x' P_z  \ \ p_y =  y' P_z
\end{equation}
$$

with $$P_z = P_0 (1+\delta)$$.  $$P_0$$ is the reference momentum.  $$\tau$$ is the time lag relative to
the ideal particle.

[1]: matlab.html "Matlab interface"
[2]: python.html "python interfqce"
