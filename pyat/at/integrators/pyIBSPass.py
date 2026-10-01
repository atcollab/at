"""Intra-beam scattering pass method, see :py:class:`.IBSElement`."""


def trackFunction(rin, elem=None):
    elem.track_turn(rin)
