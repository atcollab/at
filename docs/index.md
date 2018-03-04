# ![Image](https://cdn.rawgit.com/atcollab/atdoc/40b8230c/art/atcollab.png) Welcome to the Accelerator Toolbox WEB Pages ![Image](https://cdn.rawgit.com/atcollab/atdoc/40b8230c/art/atcollab.png)

Automatic HTML documentation is generated from mfile using the [m2html package](https://www.artefact.tk/software/matlab/m2html/) by Guillaume Flandin.

**[Online documentation of all functions](https://cdn.rawgit.com/atcollab/atdoc/aa9b9f58/doc_html/index.html)**


Introduction
============

Accelerator Toolbox (AT) is a collection of tools to model storage rings and beam transport lines in Matlab.

With AT it is possible to:

- **create and manipulate accelerator lattice elements**, which are stored in Matlab structures
- **track particles through the lattice**, selecting the appropriate integrator to represent the physics
- **compute accelerator parameters and beam properties**, generating new scripts or taking advantage of the existing ones

Examples of things AT enables are: create a lattice structure that links to the integrators, and then to track through the lattice 
including non-linear elements, analyze the non-linear motion, change settings to affect parameters, compute equilibrium beam sizes 
with radiation, and so on.

History
-------

The core of AT was developed at SLAC by Andrei Terebilo and a webpage is hosted `here <http://www.slac.stanford.edu/grp/ssrl/spear/at/>`.

The AT code in this repository evolved from version 1.3 of the original code and is now maintained by the 'atcollab' collaboration, 
involving people from different research institutes. The latest release can be found `here <https://github.com/atcollab/at/releases>`_.

Other main developpers and maintainers: 
Gregory Portmann, Laurent S. Nadolski, Eugene Tan, Xiabio Huang, C. Steier (ALS)

From version 2.0: 
Laurent Farvacque (ESRF), Simone Liuzzo (ESRF), Nicola Carminani (ESRF), Boaz Nash (ESRF), 
G. Campogiani (LNF),Laurent S. Nadolski (SOLEIL), Aurelien Bence (SOLEIL), Peace Chang (TLS), 
M. Munoz (ALBA), Z. Marti ALBA), Will Rogers (Diamond)


