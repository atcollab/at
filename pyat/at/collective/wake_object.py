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


class WakeType(Enum):
    """Enum class for wake type"""
    FILE = 1        #: Import from file
    TABLE = 2       #: Provide vectors
    RESONATOR = 3   #: Analytical resonator
    RESWALL = 4     #: Analytical resistive wall


class WakeComponent(Enum):
    """Enum class for wake component"""
    DX = 1          #: Dipole X
    DY = 2          #: Dipole Y
    QX = 3          #: Quadrupole X
    QY = 4          #: Quadrupole Y
    Z = 5           #: Longitudinal


# noinspection PyPep8Naming
class Wake(object):
    # noinspection PyUnresolvedReferences
    """Class to generate a wake object

    The wake object is defined by its srange, specified
    at initialization, and one or several components corresponding
    to transverse dipoles or quadrupoles, or longitudinal wakes.

    The ``srange`` is common to all components and cannot be changed
    once initialized, all added component are resampled to the
    ``srange``.

    Parameters:
        srange:         vector of s position where to sample the wake

    Example:
        >>> # Build the range
        >>> srange = Wake.build_srange(-0.36, 0.36, 1.0e-5, 1.0e-2,
        ... ring.circumference, ring.circumference)
        >>> # Create the Wake object
        >>> wake = Wake(srange)
        >>> # Build a resonator and add it to the Wake
        >>> freq = 10e9
        >>> Q = 1
        >>> Rs = 1e4
        >>> wake.add(WakeType.RESONATOR, WakeComponent.Z, freq, Q, Rs,
        ... ring.beta)

        See :py:meth:`Wake.long_resonator` for a simple way of creatong the
        same :py:class:`Wake`

    Components may be retrieved using :code:`wake.DX` for example
    """
    def __init__(self, srange):
        assert len(srange) == len(numpy.unique(srange)), \
            "srange must contain unique values"            
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
        """Dipole X component"""
        return self.components[WakeComponent.DX]

    @property
    def DY(self):
        """Dipole Y component"""
        return self.components[WakeComponent.DY]

    @property
    def QX(self):
        """Quadrupole X component"""
        return self.components[WakeComponent.QX]

    @property
    def QY(self):
        """Quadrupole Y component"""
        return self.components[WakeComponent.QY]

    @property
    def Z(self):
        """Longitudinal component"""
        return self.components[WakeComponent.Z]

    def add(self, wtype: WakeType, wcomp: WakeComponent, *args, **kwargs):
        """Add a component to a :py:class:`.Wake`

        Parameters:
            wtype:      Wake type
            wcomp:      Wake component
            *args:      parameters for the selected component
        """
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
        
        if numpy.amin(self._srange) < 0:
            raise ValueError("""
                            Provided srange has negative values.
                            This is not allowed for the resonator
                            wake function. Please correct.
                            """)

        # It is needed to have a point at 0 and 1e-24 to sample properly
        # the discontinuity in the longitudinal plane

        self._srange = numpy.unique(numpy.concatenate(([0.0, 1e-24],
                                                        self._srange)))

        if wcomp is WakeComponent.Z:
            return long_resonator_wf(self._srange, frequency,
                                     qfactor, rshunt, beta)
        elif isinstance(wcomp, WakeComponent):
            return transverse_resonator_wf(self._srange, frequency,
                                           qfactor, rshunt, yokoya_factor,
                                           beta)
        else:
            raise AtError('Invalid WakeComponent: {}'.format(wcomp))

    def _reswall(self, wcomp, length, rvac, conduct, beta, yokoya_factor=1):
    
        if numpy.amin(self._srange) < 0:
            raise ValueError("""
                            Provided srange has negative values.
                            This is not allowed for the resistive wall
                            wake function. Please correct.
                            """)

        # It is needed to have a point at 0 and 1e-24 to sample properly
        # the discontinuity in the longitudinal plane (if combining with
        # longres)
        self._srange = numpy.unique(numpy.concatenate(([0.0, 1e-24],
                                                       self._srange)))

        if wcomp is WakeComponent.Z:
            raise AtError('Resitive wall not available '
                          'for WakeComponent: {}'.format(wcomp))
        elif isinstance(wcomp, WakeComponent):
            return transverse_reswall_wf(self._srange, yokoya_factor,
                                         length, rvac, conduct, beta)
        else:
            raise AtError('Invalid WakeComponent: {}'.format(wcomp))

    @staticmethod
    def resonator(srange, wakecomp,
                  frequency, qfactor, rshunt,
                  beta: float, yokoya_factor=1, nelems: int = 1):
        r"""Factory method to build a resonator wake object

        Parameters:
            srange:         vector of s position where to sample the wake
            wakecomp:       Wake component
            frequency:      Resonator frequency [Hz]
            qfactor:        Q factor
            rshunt:         Shunt impedance, [:math:`\Omega`] for longitudinal,
              [:math:`\Omega/m`] for transverse
            beta:           Relativistic :math:`\beta`
            yokoya_factor:  Yokoya factor
            nelems:         Number of resonators

        Several resonators can be built in one step by setting ``nelems``> 1
        and giving array parameters
        """
        wake = Wake(srange)
        try:
            wakecomp = numpy.broadcast_to(wakecomp, (nelems, ))
            frequency = numpy.broadcast_to(frequency, (nelems, ))
            qfactor = numpy.broadcast_to(qfactor, (nelems, ))
            rshunt = numpy.broadcast_to(rshunt, (nelems, ))
            yokoya_factor = numpy.broadcast_to(yokoya_factor, (nelems, ))
        except ValueError:
            raise AtError('Wake object inputs should be either scalars '
                          'or with shape (len(wakecomp), )')
        for wc, fr, qf, rs, yk in zip(wakecomp, frequency, qfactor, rshunt,
                                      yokoya_factor):
            wake.add(WakeType.RESONATOR, wc, fr, qf, rs, beta,
                     yokoya_factor=yk)
        return wake

    @staticmethod
    def long_resonator(srange, frequency, qfactor, rshunt, beta, nelems=1):
        # noinspection PyUnresolvedReferences
        r"""Factory method to build a longitudinal resonator wake object

        Parameters:
            srange:         vector of s position where to sample the wake
            frequency:      Resonator frequency [Hz]
            qfactor:        Q factor
            rshunt:         Shunt impedance, [:math:`\Omega`] for longitudinal,
              [:math:`\Omega/m`] for transverse
            beta:           Relativistic :math:`\beta`
            nelems:         Number of resonators

        Several resonators can be built in one step by setting ``nelems``> 1
        and giving array parameters

        Example:
            >>> # Build the range
            >>> srange = Wake.build_srange(-0.36, 0.36, 1.0e-5, 1.0e-2,
            ... ring.circumference, ring.circumference)
            >>> # Create the resonator
            >>> freq = 10e9
            >>> Q = 1
            >>> Rs = 1e4
            >>> wake = Wake.long_resonator(srange, freq, Q, Rs, ring.beta)
        """ 
        return Wake.resonator(srange, WakeComponent.Z, frequency, qfactor,
                              rshunt, beta, nelems=nelems)

    @staticmethod
    def resistive_wall(srange, wakecomp: WakeComponent,
                       length, rvac, conduct, beta: float,
                       yokoya_factor=1, nelems: int = 1):
        r"""Factory method to build a resistive wall wake object

        Parameters:
            srange:         vector of s position where to sample the wake
            wakecomp:       Wake component
            length:
            rvac:
            conduct:
            beta:           Relativistic :math:`\beta`
            yokoya_factor:  Yokoya factor
            nelems:
        """
        wake = Wake(srange)
        try:
            wakecomp = numpy.broadcast_to(wakecomp, (nelems, ))
            length = numpy.broadcast_to(length, (nelems, ))
            rvac = numpy.broadcast_to(rvac, (nelems, ))
            conduct = numpy.broadcast_to(conduct, (nelems, ))
            yokoya_factor = numpy.broadcast_to(yokoya_factor, (nelems, ))
        except ValueError:
            raise AtError('Wake object inputs should be either scalars '
                          'or with shape (len(wakecomp), )')
        for wc, le, rv, co, yk in zip(wakecomp, length, rvac,
                                      conduct, yokoya_factor):
            wake.add(WakeType.RESWALL, wc, le, rv, co, beta, yokoya_factor=yk)
        return wake

    @staticmethod
    def build_srange(start: float, bunch_ext: float, short_step: float,
                     long_step: float,
                     bunch_interval: float, totallength: float):
        """Function to build the wake table s column.
        This is not the slicing but the look-up table,
        however it generates data where bunches are located
        to avoid using too much memory to store the table.

        Parameters:
            start:          starting s-coordinate of the table
                            (can be negative for wake potential)
            bunch_ext:      maximum bunch extension, function
                            generates data at +/- bunch_ext
                            around the bucket center
            short_step:     step size for the short range wake table
            long_step:      step size for the long range wake table
            bunch_interval: minimum bunch interval. Data will be generated
                            for each ``bunch_interval`` step
            totallength:    total length of the wake table, has to contain
                            the full bunch extension

        Returns:
            srange:         vector of s position where to sample the wake
        """
        srange = numpy.arange(start, bunch_ext, short_step)
        rangel = numpy.arange(-bunch_ext, bunch_ext, long_step)
        nbunch = int((totallength-bunch_ext)/bunch_interval)
        for i in range(nbunch):
            srange = numpy.concatenate((srange, rangel+bunch_interval*(i+1)))
        return numpy.unique(srange)
