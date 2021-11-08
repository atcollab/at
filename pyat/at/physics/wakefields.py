"""
Wake field and wake element generation
"""
import at
import numpy
import warnings
from enum import Enum
from scipy.constants import physical_constants
from scipy.constants import c as clight
from scipy.constants import e as qe
from scipy.interpolate import interp1d
from at.lattice import AtWarning, AtError
partmass = physical_constants['electron mass energy'
                              ' equivalent in MeV'][0]*1.0e6


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
    RESWALL = 4 # Analytical resistive wall


class WakeComponent(Enum):
    """Enum class for wake component"""
    DX = 1  # Dipole X
    DY = 2  # Dipole Y
    QX = 3  # Quadrupole X
    QY = 4  # Quadrupole Y
    Z = 5   # Longitudinal


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
    (analytical resonator)

    Components are retrieved with Wake.DX for example
    """
    def __init__(self, srange):
        self._srange = srange
        self.components = {WakeComponent.DX: None,
                           WakeComponent.DY: None,
                           WakeComponent.QX: None,
                           WakeComponent.QY: None,
                           WakeComponent.Z: None}

    # noinspection PyPep8Naming
    @property
    def DX(self):
        return self.components[WakeComponent.DX]

    # noinspection PyPep8Naming
    @DX.setter
    def DX(self, s, w):
        self.components[WakeComponent.DX] = self.resample(s, w)

    # noinspection PyPep8Naming
    @property
    def DY(self):
        return self.components[WakeComponent.DY]

    # noinspection PyPep8Naming
    @DY.setter
    def DY(self, s, w):
        self.components[WakeComponent.DY] = self.resample(s, w)

    # noinspection PyPep8Naming
    @property
    def QX(self):
        return self.components[WakeComponent.QX]

    # noinspection PyPep8Naming
    @QX.setter
    def QX(self, s, w):
        self.components[WakeComponent.QX] = self.resample(s, w)

    # noinspection PyPep8Naming
    @property
    def QY(self):
        return self.components[WakeComponent.QY]

    # noinspection PyPep8Naming
    @QY.setter
    def QY(self, s, w):
        self.components[WakeComponent.QY] = self.resample(s, w)

    # noinspection PyPep8Naming
    @property
    def Z(self):
        return self.components[WakeComponent.Z]

    # noinspection PyPep8Naming
    @Z.setter
    def Z(self, s, w):
        self.components[WakeComponent.Z] = self.resample(s, w)

    def get_srange(self):
        return self._srange

    def add(self, wtype, wcomp, *args, **kwargs):
        if wtype is WakeType.FILE:
            w = self.readwakefile(*args, **kwargs)
        elif wtype is WakeType.TABLE:
            w = self.resample(*args)
        elif wtype is WakeType.RESONATOR:
            w = self.resonator(wcomp, *args, **kwargs)
        else:
            raise AtError('Invalid WakeType: {}'.format(wtype))
        if self.components[wcomp] is None:
            self.components[wcomp] = w
        else:
            self.components[wcomp] += w

    def resample(self, s, w):
        if self._srange[0] < s[0] or self._srange[-1] < s[-1]:
            warnings.warn(AtWarning('Input wake is smaller '
                                    'than desired Wake() range. '
                                    'Filling with zeros.\n'))
        fint = interp1d(s, w, bounds_error=False, fill_value=0)
        wint = fint(self._srange)
        return wint

    def convolve_wakefun(self, w, sigt):        
        min_step = numpy.diff(self._srange)
        t_out = numpy.arange(self._srange[0], self._srange[-1], min_step)
        sdiff = t_out[-1]-t_out[0]
        npoints = len(t_out)
        nt = npoints+npoints-1;
        func = interp1d(self._srange, w, bounds_error=False, fill_value = 0)
        wout = func(t_out)
        wout = numpy.append(wout, numpy.zeros(nt-len(wout)))
        fftr = numpy.fft.fft(wout)
        f = numpy.fft.fftshift(np.linspace(-(npoints-1)/sdiff,(npoints-1)/sdiff,nt))
        fftl = numpy.exp(-(f*2*np.pi*sigt)**2/2)
        wout = numpy.fft.ifft(fftr*fftl)
        wout = numpy.roll(wout,int(npoints/2))
        t_out = numpy.linspace(t_out[0],t_out[-1],nt)
        func = interp1d(t_out, wout, bounds_error=False, fill_value = 0)
        wout = func(self._srange)
    return wout

    def readwakefile(self, filename, scol=0, wcol=1, sfact=1, wfact=1,
                     delimiter=None, skiprows=0):
        s, w = numpy.loadtxt(filename, delimiter=delimiter, unpack=True,
                             usecols=(scol, wcol), skiprows=skiprows)
        s *= sfact
        w *= wfact
        return self.resample(s, w)

    def resonator(self, wcomp, frequency, qfactor, rshunt, beta,
                  yokoya_factor=1):
        if wcomp is WakeComponent.Z:
            return self.wakefunc_long_resonator(frequency, qfactor,
                                                rshunt, beta)
        elif wcomp is WakeComponent.DX:
            return self.wakefunc_transverse_resonator(frequency, qfactor,
                                                      rshunt, yokoya_factor,
                                                      beta)
        elif wcomp is WakeComponent.DY:
            return self.wakefunc_transverse_resonator(frequency, qfactor,
                                                      rshunt, yokoya_factor,
                                                      beta)
        elif wcomp is WakeComponent.QX:
            return self.wakefunc_transverse_resonator(frequency, qfactor,
                                                      rshunt, yokoya_factor,
                                                      beta)
        elif wcomp is WakeComponent.QY:
            return self.wakefunc_transverse_resonator(frequency, qfactor,
                                                      rshunt, yokoya_factor,
                                                      beta)
        else:
            raise AtError('Invalid WakeComponent: {}'.format(wcomp))

    def reswall(self, wcomp, length, rvac, conduct, beta, yokoya_factor=1):
        if wcomp is WakeComponent.Z:
            raise AtError('Resitive wall not available '
                          'for WakeComponent: {}'.format(wcomp))
        elif wcomp is WakeComponent.DX:
            return self.wakefunc_reswall(yokoya_factor, length, 
                                         rvac, conduct, beta)
        elif wcomp is WakeComponent.DY:
            return self.wakefunc_reswall(yokoya_factor, length,
                                         rvac, conduct, beta)
        elif wcomp is WakeComponent.QX:
            raise AtError('Resitive wall not available '
                          'for WakeComponent: {}'.format(wcomp))
        elif wcomp is WakeComponent.QY:
            raise AtError('Resitive wall not available '
                          'for WakeComponent: {}'.format(wcomp))
        else:
            raise AtError('Invalid WakeComponent: {}'.format(wcomp))

    def wakefunc_long_resonator(self, frequency, qfactor, rshunt, beta):
        """Define the wake function (longitudinal) of a resonator
        with the given parameters according to Alex Chao's resonator
        model (Eq. 2.82) and definitions of the resonator in HEADTAIL.
        """
        omega = 2 * numpy.pi * frequency
        alpha = omega / (2 * qfactor)
        omegabar = numpy.sqrt(numpy.abs(omega**2 - alpha**2))
        dt = -self._srange/(beta * clight)
        if qfactor > 0.5:
            wake = (-(numpy.sign(dt) - 1) * rshunt * alpha *
                    numpy.exp(alpha * dt) * (numpy.cos(omegabar * dt) +
                    alpha / omegabar * numpy.sin(omegabar*dt)))
        elif qfactor == 0.5:
            wake = (-(numpy.sign(dt) - 1) * rshunt * alpha *
                    numpy.exp(alpha * dt) * (1. + alpha * dt))
        elif qfactor < 0.5:
            wake = (-(numpy.sign(dt) - 1) * rshunt * alpha *
                    numpy.exp(alpha * dt) * (numpy.cosh(omegabar * dt) +
                    alpha / omegabar * numpy.sinh(omegabar * dt)))
        return wake

    def wakefunc_transverse_resonator(self, frequency, qfactor, rshunt,
                                      yokoya_factor, beta):
        """Define the wake function (transverse) of a resonator
        with the given parameters according to Alex Chao's
        resonator model (Eq. 2.82) and definitions of the resonator
        in HEADTAIL.
        """
        omega = 2 * numpy.pi * self.frequency
        alpha = omega / (2 * qfactor)
        omegabar = numpy.sqrt(numpy.abs(omega**2 - alpha**2))
        dt = -self._srange/(beta * clight)
        if qfactor > 0.5:
            wake = (yokoya_factor * rshunt * omega**2 / (qfactor *
                    omegabar) * numpy.exp(alpha*dt) * numpy.sin(omegabar*dt))
        elif qfactor == 0.5:
            wake = (yokoya_factor * rshunt * omega**2 / qfactor *
                    numpy.exp(alpha * dt) * dt)
        else:
            wake = (yokoya_factor * rshunt * omega**2 / (qfactor *
                    omegabar) * numpy.exp(alpha*dt) * numpy.sinh(omegabar*dt))
        return wake

    def wakefunc_reswall(self, yokoya_factor, length, rvac, conduct, beta):
        """Define the wake function (transverse) of a resistive wall with the given
        parameters according to Alex Chao's RW model (2.53) and definitions used in
        HEADTAIL
        """
        z0 = 119.9169832 * np.pi
        dt = -self._srange/(beta * clight)
        wake = yokoya_factor * (numpy.sign(dt) - 1) / 2. *
               beta * length / numpy.pi / rvac **3*
               numpy.sqrt(-clight * z0 /conduct / numpy.pi / dt))
        return wake


class WakeElement(at.Element):
    """Class to generate an AT wake element using the
    passmethod WakeFieldPass
    args:  family name, ring, wake object
    kwargs: Intensity  (default=0) bunch intensity
            Passmethod (default=WakeFieldPass)
            Nslice     (default=101) number of slices
                       per bunch
            Nturns     (default=1) number of turn for 
                       the wake field
            ZCuts      (default=None)limits for fixed 
                       slicing, default is adaptative
            NormFact   (default=[1,1,1]) normalization
                       for the 3 planes, to account for
                       beta function at the observation
                       point for example
    """
    def __init__(self, family_name, ring, wake, **kwargs):         
        kwargs.setdefault('PassMethod', 'WakeFieldPass')       
        self.Intensity = kwargs.pop('Intensity', 0.0)
        self.Nslice = kwargs.pop('Nslice', 101)
        self.NormFact = kwargs.pop('NormFact',
                                   numpy.ones(3,order='F'))
        self.Wakefact = self.get_wakefact(ring)
        self.int2curr = self.get_int2curr(ring)
        self.WakeT = wake.get_srange()
        self.Nelem = len(self.WakeT)
        zcuts = kwargs.pop('ZCuts',None)
        if zcuts is not None:
            self.ZCuts=zcuts
        if wake.Z is not None:
            self.WakeZ = wake.Z
        if wake.DX is not None:
            self.WakeDX = wake.DX
        if wake.DY is not None:
            self.WakeDY = wake.DY
        if wake.QX is not None:
            self.WakeQX = wake.QX
        if wake.QY is not None:
            self.WakeQY = wake.QY
        self.Nturns = kwargs.pop('Nturns', 1)
        self.Circumference = ring.circumference
        self.TurnHistory = numpy.zeros((self.Nturns*self.Nslice,4),
                                       order='F')
        super(WakeElement, self).__init__(family_name, **kwargs)      

    def get_wakefact(self, ring):
        betrel = numpy.sqrt(1.0-(partmass/ring.energy)**2)
        return -qe/(ring.energy*betrel**2)

    def get_int2curr(self, ring):
        betrel = numpy.sqrt(1.0-(partmass/ring.energy)**2)
        return clight*betrel*qe/ring.circumference

    def clear_history(self):
        self.TurnHistory = numpy.zeros((self.Nturns*self.Nslice,4),
                                       order='F')

    # noinspection PyPep8Naming
    @property
    def Current(self):
        return self.Intensity*self.int2curr

    # noinspection PyPep8Naming
    @Current.setter
    def Current(self, current):
        self.Intensity = current/self.int2curr

    def __repr__(self):
        """Simplified __repr__ to avoid errors due to arguments 
        not defined as attributes
        """
        attrs = dict(self.items())
        return '{0}({1})'.format(self.__class__.__name__, attrs)


class LongResonatorElement(WakeElement):
    """Class to generate a longitudinal resonator, inherits from Wake()
       additonal argument are frequency, qfactor, rshunt
    """
    def __init__(self, family_name, ring, srange, frequency, qfactor,
                 rshunt, **kwargs):
        self.Frequency = frequency
        self.QFactor = qfactor
        self.RShunt = rshunt
        beta = numpy.sqrt(1.0-(partmass/ring.energy)**2)
        wake = Wake(srange)
        wake.add(WakeType.RESONATOR, WakeComponent.Z, frequency,
                 qfactor, rshunt, beta)
        super(LongResonatorElement, self).__init__(family_name, ring=ring,
                                                   wake=wake, **kwargs)
