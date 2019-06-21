"""Plot a lattice synoptic"""
import numpy
# noinspection PyPackageRequirements
import matplotlib.pyplot as plt
# noinspection PyPackageRequirements
from matplotlib.patches import Polygon
# noinspection PyPackageRequirements
from matplotlib.collections import PatchCollection
from at.lattice import elements as elts

__all__ = ['plot_synopt']

# Default properties for element representation
DIPOLE = dict(label='Dipoles', facecolor=(0.5, 0.5, 1.0))
QUADRUPOLE = dict(label='Quadrupoles', facecolor=(1.0, 0.5, 0.5))
SEXTUPOLE = dict(label='Sextupoles', facecolor=(0.5, 1.0, 0.5))
MULTIPOLE = dict(label='Multipoles', facecolor=(0.25, 0.75, 0.25))
MONITOR = dict(label='Monitors', linestyle=None, marker=10, color='k')


# noinspection PyDefaultArgument
def plot_synopt(ring, axes=None, dipole={}, quadrupole={}, sextupole={},
                multipole={}, monitor={}):
    """Plot a synoptic of a lattice

    PARAMETERS
        ring            Lattice object

    KEYWORDS
        s_range=None    plot range, defaults to the full ring
        axes=None       axes for plotting the synoptic. If None, a new
                        figure will be created. Otherwise, a new axes object
                        sharing the same x-axis as the given one is created.
        dipole={}       Dictionary of properties overloading the default
                        properties. If None, dipoles will not be shown.
        quadrupole={}   Same definition as for dipole
        sextupole={}    Same definition as for dipole
        multipole={}    Same definition as for dipole
        monitor={}      Same definition as for dipole

    RETURN
        synopt_axes     Synoptic axes
     """

    class Dipole(Polygon):
        xx = numpy.array([0, 0, 1, 1], dtype=float)
        yy = numpy.array([0, 1, 1, 0], dtype=float)

        def __init__(self, s, length, **kwargs):
            xy = numpy.stack((self.xx * length + s, self.yy), axis=1)
            super(Dipole, self).__init__(xy, closed=False, **kwargs)

    class Quadrupole(Polygon):
        xx = numpy.array([0, 0, 0.5, 1, 1])
        yy = {True: numpy.array([0, 1, 1.4, 1, 0]),
              False: numpy.array([0, 1, 0.6, 1, 0])}

        def __init__(self, s, length, foc, **kwargs):
            xy = numpy.stack((self.xx * length + s, self.yy[foc]), axis=1)
            super(Quadrupole, self).__init__(xy, closed=False, **kwargs)

    class Sextupole(Polygon):
        xx = numpy.array([0, 0, 0.33, 0.66, 1, 1])
        yy = {True: numpy.array([0, 0.8, 1, 1, 0.8, 0]),
              False: numpy.array([0, 0.8, 0.6, 0.6, 0.8, 0])}

        def __init__(self, s, length, foc, **kwargs):
            xy = numpy.stack((self.xx * length + s, self.yy[foc]), axis=1)
            super(Sextupole, self).__init__(xy, closed=False, **kwargs)

    class Multipole(Polygon):
        xx = numpy.array([0, 0, 1, 1], dtype=float)
        yy = numpy.array([0, 0.8, 0.8, 0])

        def __init__(self, s, length, **kwargs):
            xy = numpy.stack((self.xx * length + s, self.yy), axis=1)
            super(Multipole, self).__init__(xy, closed=False, **kwargs)

    class Monitor(Polygon):
        xx = numpy.array([0.0, 0.0])
        yy = numpy.array([0.0, 1.2])

        def __init__(self, s, **kwargs):
            xy = numpy.stack((self.xx + s, self.yy), axis=1)
            super(Monitor, self).__init__(xy, closed=False, **kwargs)

    def ismultipole(elem):
        return isinstance(elem, elts.Multipole) and not isinstance(elem, (
            elts.Dipole, elts.Quadrupole, elts.Sextupole))

    if axes is None:
        fig = plt.figure()
        axsyn = fig.add_subplot(111, xlim=ring.s_range)
    else:
        axsyn = axes.twinx()
    axsyn.set_axis_off()         # Set axis invisible
    axsyn.set_ylim((0.0, 20.0))  # Initial scaling of elements
    axsyn.set_zorder(-0.2)       # Put synoptic in the background

    s_pos = ring.get_s_pos(range(len(ring)))

    if dipole is not None:
        props = DIPOLE.copy()
        props.update(dipole)
        dipoles = PatchCollection(
            (Dipole(s, el.Length) for s, el in zip(s_pos, ring)
             if isinstance(el, elts.Dipole)), **props)
        axsyn.add_collection(dipoles)

    if quadrupole is not None:
        props = QUADRUPOLE.copy()
        props.update(quadrupole)
        quadrupoles = PatchCollection(
            (Quadrupole(s, el.Length, el.PolynomB[1] >= 0.0)
             for s, el in zip(s_pos, ring)
             if isinstance(el, elts.Quadrupole)), **props)
        axsyn.add_collection(quadrupoles)

    if sextupole is not None:
        props = SEXTUPOLE.copy()
        props.update(sextupole)
        sextupoles = PatchCollection(
            (Sextupole(s, el.Length, el.PolynomB[2] >= 0.0)
             for s, el in zip(s_pos, ring)
             if isinstance(el, elts.Sextupole)), **props)
        axsyn.add_collection(sextupoles)

    if multipole is not None:
        props = MULTIPOLE.copy()
        props.update(multipole)
        multipoles = PatchCollection(
            (Multipole(s, el.Length) for s, el in zip(s_pos, ring)
             if ismultipole(el)), **props)
        axsyn.add_collection(multipoles)

    if monitor is not None:
        props = MONITOR.copy()
        props.update(monitor)
        s = s_pos[[isinstance(el, elts.Monitor) for el in ring]]
        y = numpy.zeros(s.shape)
        # noinspection PyUnusedLocal
        monitors = axsyn.plot(s, y, **props)

    return axsyn
