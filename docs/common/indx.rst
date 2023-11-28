AT basics
=========

The Accelerator Toolbox is based on a fast and modular tracking engine.

Coordinate system
-----------------
AT runs in 6 dimensions using the following coordinate system:

.. math::

   \vec Z =
   \begin{pmatrix}
   x \\ p_x \\
   y \\ p_y \\
   \delta \\ \beta c\tau
   \end{pmatrix}

AT works with relative path lengths: the 6\ :sup:`th` coordinate :math:`\beta c\tau`
represents the path lengthening with respect to the reference particle. :math:`\tau`
is the delay with respect to the reference particle: the particle is late for
:math:`\tau > 0`.

The reference particle is conventionally a particle travelling on the reference orbit
at constant nominal energy and crossing the RF cavities at their 0 voltage.

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
