AT basics
=========

The Accelerator Toolbox is based on a fast and modular tracking engine.

Coordinate system
-----------------
The AT coordinate system is based on a design trajectory along which the magnetic
elements are aligned, and a reference particle:

.. image:: /images/coord_system.*
   :alt: coordinate system
   :align: center

By convention, the reference particle travels along the design trajectory at constant
nominal energy and defines the phase origin of the RF voltage.

The 6-d coordinate vector :math:`\vec Z` is:

.. math::

   \vec Z =
   \begin{pmatrix}
   x \\ p_x/p_0 \\
   y \\ p_y/p_0 \\
   \delta \\ \beta c\tau
   \end{pmatrix}

:math:`p_0` is the nominal momentum,

:math:`\delta=\dfrac{p-p_0}{p_0}` is the relative momentum deviation.

AT works with relative path lengths: the 6\ :sup:`th` coordinate :math:`\beta c\tau`
represents the path lengthening with respect to the reference particle. :math:`\tau`
is the delay with respect to the reference particle: the particle is late for
:math:`\tau > 0`.

Tracking
--------
All the AT processing is based on tracking. For instance, linear optics properties
are deduced from a one-turn transfer matrix obtained by differentiating the tracking
output on a small dimension grid centered on the closed orbit.

The tracking through each AT element is performed by an integrator, called its
"passmethod". There is no correlation between an element class and its associated
passmethod. The passmethod must be explicitly specified by the ``PassMethod``
attribute of the element.

The tracking engine is written in C for the best performance. Each passmethod
is dynamically loaded on-demand when encountering the first element requiring it.

See :doc:`passmethods` for a description of the available passmethods,

.. toctree::
   :hidden:

   Passmethods <passmethods>
   About AT <about>
