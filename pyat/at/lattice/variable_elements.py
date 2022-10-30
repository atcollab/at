import numpy
from .elements import Element, _array
from .utils import AtError
from enum import IntEnum
from datetime import datetime


class ACMode(IntEnum):
    SINE = 0
    WHITENOISE = 1
    ARBITRARY = 2


class VariableMultipole(Element):
    _BUILD_ATTRIBUTES = Element._BUILD_ATTRIBUTES
    _conversions = dict(Element._conversions,
                        Mode=int,
                        AmplitudeA=_array, AmplitudeB=_array,
                        FrequencyA=float, FrequencyB=float,
                        PhaseA=float, PhaseB=float,
                        Seed=int, NSamplesA=int, NSamplesB=int,
                        FuncA=_array, FuncB=_array, Ramps=_array)

    def __init__(self, family_name, AmplitudeA=None, AmplitudeB=None,
                 mode=ACMode.SINE, MaxOrder=0, **kwargs):
        kwargs.setdefault('PassMethod', 'VariableThinMPolePass')
        self.MaxOrder = MaxOrder
        self.Mode = int(mode)
        if AmplitudeA is None and AmplitudeB is None:
            raise AtError('Please provide at least one amplitude for A or B')
        AmplitudeB = self._set_params(AmplitudeB, mode, 'B', **kwargs)
        AmplitudeA = self._set_params(AmplitudeA, mode, 'A', **kwargs)
        self._setmaxorder(AmplitudeA, AmplitudeB)
        if mode == ACMode.WHITENOISE:
            self.Seed = kwargs.pop('Seed', datetime.now().timestamp())
        self.PolynomA = numpy.zeros(self.MaxOrder+1)
        self.PolynomB = numpy.zeros(self.MaxOrder+1)
        ramps = kwargs.pop('Ramps', None)
        if ramps is not None:
            self.Ramps = ramps
        super(VariableMultipole, self).__init__(family_name, **kwargs)

    def _setmaxorder(self, ampa, ampb):
        mxa, mxb = 0, 0
        if ampa is not None:
            mxa = numpy.max(numpy.nonzero(ampa))
        if ampb is not None:
            mxb = numpy.max(numpy.nonzero(ampb))
        self.MaxOrder = max(mxa, mxb)
        if ampa is not None:
            delta = self.MaxOrder - len(ampa)
            if delta > 0:
                ampa = numpy.pad(ampa, (0, delta))
            setattr(self, 'AmplitudeA', ampa)
        if ampb is not None:
            delta = self.MaxOrder + 1 - len(ampb)
            if delta > 0:
                ampb = numpy.pad(ampb, (0, delta))
            setattr(self, 'AmplitudeB', ampb)

    def _set_params(self, amplitude, mode, ab, **kwargs):
        if amplitude is not None:
            if numpy.isscalar(amplitude):
                amp = numpy.zeros(self.MaxOrder)
                amplitude = numpy.append(amp, amplitude)
            if mode == ACMode.SINE:
                self._set_sine(ab, **kwargs)
            if mode == ACMode.ARBITRARY:
                self._set_arb(ab, **kwargs)
        return amplitude

    def _set_sine(self, ab, **kwargs):
        frequency = kwargs.pop('Frequency' + ab, None)
        phase = kwargs.pop('Phase' + ab, 0)
        assert frequency is not None, \
            'Please provide a value for Frequency' + ab
        setattr(self, 'Frequency' + ab, frequency)
        setattr(self, 'Phase' + ab, phase)

    def _set_arb(self, ab, **kwargs):
        func = kwargs.pop('Func' + ab, None)
        nsamp = len(func)
        assert func is not None, \
            'Please provide a value for Func' + ab
        setattr(self, 'Func' + ab, func)
        setattr(self, 'NSamples' + ab, nsamp)
