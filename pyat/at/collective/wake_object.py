"""
Wake object creation
"""
import numpy
import warnings
from enum import Enum
from scipy.interpolate import interp1d
from ..lattice import AtWarning, AtError
from .wake_functions import long_resonator_wf, transverse_resonator_wf
from .wake_functions import transverse_reswall_wf


def build_srange(start, bunch_ext, short_step, long_step,
                 bunch_interval, totallength):
    """Function to build the wake table s column.
    This is not the slicing but the look-up table,
    however it generates data where bunches are located
    to avoid using too much memory to store the table.

    PARAMETERS
        start           starting s-coordinate of the table
                        (can be negative for wake potential)
        bunch_ext       maximum bunch extension, function
                        generates data at +/- bunch_ext
                        around the bucket center
        short_step      step size for the short range wake table
        long_step       step size for the long range wake table
        bunch_interval  minimum bunch interval data will be generate
                        for each bunch_inteval step
        totallength     total length of the wake table, has to contain
                        the full bunch extension

    OUTPUT
        srange          vector of s position where to sample the wake
    """
    srange = numpy.arange(start, bunch_ext, short_step)
    rangel = numpy.arange(-bunch_ext, bunch_ext, long_step)
    nbunch = int((totallength-bunch_ext)/bunch_interval)
    for i in range(nbunch):
        srange = numpy.concatenate((srange, rangel+bunch_interval*(i+1)))
    return numpy.unique(srange)


class WakeType(Enum):
    """Enum class for wake type"""
    FILE = 1  # Import from file
    TABLE = 2  # Provide vectors
    RESONATOR = 3  # Analytical resonator
    RESWALL = 4  # Analytical resistive wall


class WakeComponent(Enum):
    """Enum class for wake component"""
    DX = 1  # Dipole X
    DY = 2  # Dipole Y
    QX = 3  # Quadrupole X
    QY = 4  # Quadrupole Y
    Z = 5   # Longitudinal


# noinspection PyPep8Naming
class Wake(object):
    """Class to generate a wake object
    The wake object is define by its srange, specified
    at initialization, and DX, DY, QY, Z corresponding
    to transverse dipoles and quadrupoles and longitudinal

    The srange is common to all components and cannot be changed
    once initialized, all added component are resampled to the
    srange

    usage:
    wake = Wake(srange)
    wake.add(WakeType,WakeComponent, *args, *kwargs)

    Component are WakeComponent.FILE (import from file),
    WakeComponent.TABLE (provide vectors), WakeComponent.RESONATOR
    (analytical resonator), WakeComponent.RESWALL (transverse RW)

    Components are retrieved with Wake.DX for example
    """
    def __init__(self, srange):
        self._srange = srange
        self.components = {WakeComponent.DX: None,
                           WakeComponent.DY: None,
                           WakeComponent.QX: None,
                           WakeComponent.QY: None,
                           WakeComponent.Z: None}

    @property
    def srange(self):
        return self._srange

    @property
    def DX(self):
        return self.components[WakeComponent.DX]

    @property
    def DY(self):
        return self.components[WakeComponent.DY]

    @property
    def QX(self):
        return self.components[WakeComponent.QX]

    @property
    def QY(self):
        return self.components[WakeComponent.QY]

    @property
    def Z(self):
        return self.components[WakeComponent.Z]

    def add(self, wtype, wcomp, *args, **kwargs):
        if wtype is WakeType.FILE:
            w = self._readwakefile(*args, **kwargs)
        elif wtype is WakeType.TABLE:
            w = self._resample(*args)
        elif wtype is WakeType.RESONATOR:
            w = self._resonator(wcomp, *args, **kwargs)
        elif wtype is WakeType.RESWALL:
            w = self._reswall(wcomp, *args, **kwargs)

        else:
            raise AtError('Invalid WakeType: {}'.format(wtype))
        if self.components[wcomp] is None:
            self.components[wcomp] = w
        else:
            self.components[wcomp] += w

    def _resample(self, s, w):
        if self._srange[0] < s[0] or self._srange[-1] < s[-1]:
            warnings.warn(AtWarning('Input wake is smaller '
                                    'than desired Wake() range. '
                                    'Filling with zeros.\n'))
        fint = interp1d(s, w, bounds_error=False, fill_value=0)
        wint = fint(self._srange)
        return wint

    def _readwakefile(self, filename, scol=0, wcol=1, sfact=1, wfact=1,
                     delimiter=None, skiprows=0):
        s, w = numpy.loadtxt(filename, delimiter=delimiter, unpack=True,
                             usecols=(scol, wcol), skiprows=skiprows)
        s *= sfact
        w *= wfact
        return self._resample(s, w)

    def _resonator(self, wcomp, frequency, qfactor, rshunt, beta,
                  yokoya_factor=1):
        if wcomp is WakeComponent.Z:
            return long_resonator_wf(self._srange, frequency,
                                  qfactor, rshunt, beta)
        elif isinstance(wcomp, WakeComponent):
            return transverse_resonator_wf(self._srange, frequency,
                                        qfactor, rshunt, yokoya_factor, beta)
        else:
            raise AtError('Invalid WakeComponent: {}'.format(wcomp))

    def _reswall(self, wcomp, length, rvac, conduct, beta, yokoya_factor=1):
        if wcomp is WakeComponent.Z:
            raise AtError('Resitive wall not available '
                          'for WakeComponent: {}'.format(wcomp))
        elif isinstance(wcomp, WakeComponent):
            return transverse_reswall_wf(self._srange, yokoya_factor,
                                      length, rvac, conduct, beta)
        else:
            raise AtError('Invalid WakeComponent: {}'.format(wcomp))

    @staticmethod
    def resonator(srange, wakecomp, frequency, qfactor, rshunt, beta,
                  yokoya_factor=1):
        """
        Method to build a resonator wake object
        """
        wake = Wake(srange)
        wakecomp = numpy.atleast_1d(wakecomp)
        try:
            frequency = numpy.broadcast_to(frequency, (len(wakecomp), ))
            qfactor = numpy.broadcast_to(qfactor, (len(wakecomp), ))
            rshunt = numpy.broadcast_to(rshunt, (len(wakecomp), ))
            yokoya_factor = numpy.broadcast_to(yokoya_factor, (len(wakecomp), )) 
        except ValueError:
            raise AtError('Wake object inputs should be either scalars '
                          'or with shape (len(wakecomp), )')
        for wc, fr, qf, rs, yk in zip(wakecomp, frequency, qfactor, rshunt, yokoya_factor):
            wake.add(WakeType.RESONATOR, wc, fr, qf, rs, beta, yokoya_factor=yk)
        return wake

    @staticmethod
    def long_resonator(srange, frequency, qfactor, rshunt, beta):
        """
        Method to build a longitudinal resonator wake object
        """
        return Wake.resonator(srange, WakeComponent.Z, frequency, qfactor, rshunt, beta)

    @staticmethod
    def resistive_wall(srange, wakecomp, length, rvac, conduct, beta, yokoya_factor=1):
        """
        Method to build a resistive wall wake object
        """
        wake = Wake(srange)
        wakecomp = numpy.atleast_1d(wakecomp)
        try:
            length = numpy.broadcast_to(length, (len(wakecomp), ))
            rvac = numpy.broadcast_to(rvac, (len(wakecomp), ))
            conduct = numpy.broadcast_to(conduct, (len(wakecomp), ))
            yokoya_factor = numpy.broadcast_to(yokoya_factor, (len(wakecomp), )) 
        except ValueError:
            raise AtError('Wake object inputs should be either scalars '
                          'or with shape (len(wakecomp), )')
        for wc, le, rv, co, yk in zip(wakecomp, length, rvac, conduct, yokoya_factor):
            wake.add(WakeType.RESWALL, wc, le, rv, co, beta, yokoya_factor=yk)
        return wake
