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
                 mode=ACMode.SINE, **kwargs):
        kwargs.setdefault('PassMethod', 'VariableThinMPolePass')
        self.MaxOrder = 1
        self.Mode = int(mode)
        if AmplitudeA is None and AmplitudeB is None:
            raise AtError('Please provide at least one amplitude vector')
        self._set_params(AmplitudeB, mode, 'B', **kwargs)
        self._set_params(AmplitudeA, mode, 'A', **kwargs)
        if mode == ACMode.WHITENOISE:
            self.Seed = kwargs.pop('Seed', datetime.now().timestamp())
        self.PolynomA = numpy.zeros(self.MaxOrder)
        self.PolynomB = numpy.zeros(self.MaxOrder)
        ramps = kwargs.pop('Ramps', None)
        if ramps is not None:
            self.Ramps = ramps
        super(VariableMultipole, self).__init__(family_name, **kwargs)

    def _set_params(self, amplitude, mode, ab, **kwargs):
        if amplitude is not None:
            setattr(self, 'Amplitude' + ab, amplitude)
            if mode == ACMode.SINE:
                self._set_sine(ab, **kwargs)
            if mode == ACMode.ARBITRARY:
                self._set_arb(ab, **kwargs)
            self.MaxOrder = max(self.MaxOrder, len(amplitude))

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


class VariableMagnet(VariableMultipole):

    def __init__(self, family_name, order=1, AmplitudeA=None,
                 AmplitudeB=None, **kwargs):
        amp = numpy.zeros(order)
        if AmplitudeA is not None:
            AmplitudeA = numpy.append(amp, AmplitudeA)
            kwargs['AmplitudeA'] = AmplitudeA
        if AmplitudeB is not None:
            AmplitudeB = numpy.append(amp, AmplitudeB)
            kwargs['AmplitudeB'] = AmplitudeB
        super(VariableMagnet, self).__init__(family_name, **kwargs)
