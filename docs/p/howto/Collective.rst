Collective
==========

Overview
--------

A collective effects element can be added to any lattice to model
impedance driven collective effects and perform multi-particle tracking.
This element can be constructed with either a user defined table, or by
using the built in functions to construct specific elements (for example
a longitudinal resonator or a transverse resistive wall). The element will
then call the **WakeFieldPass** PassMethod. 

The wake field is applied by first uniformly slicing the full region occupied by the 
particles in the ct co-ordinate. Each particle is then attributed to a
given slice, which is represented by a weight. The mean position in x,y,ct 
is computed for each slice. The kick from one slice to the next (or self kick) can then be computed by taking into account
the differences in offsets of each slice. This kick is then applied to each particle 
within the slice. 

To take into account multi-turn wakes, the wake element has a **TurnHistory** buffer.
Each turn, the mean x,y,ct and weight of each slice is recorded. After each turn, the 
array of ct values is increased by one circumference (to take into account the decay 
between turns). When the kick is computed, the full history of turns is used. 

The package is organised as follows:
**at.collective.wake_functions** contains the analytic wake functions that are called  



.. math:: \\begin{equation} W(\tau) = \left\{ \begin{array}{lr} \alpha R_{s} \text{for \tau=0} \\ 2\alpha R_{s}e^{-\alpha \tau} [cos(\omega_{bar}\tau) - \frac{\alpha}{\omega_{bar}}sin(\omega_{bar}\tau] \text{for \tau > 0} \\ \end{array} \right. \\end{equation} 


