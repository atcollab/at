"""Pytohn implementation of ChangeRefPass."""

def trackFunction(rin, elem=None):
    scaling = elem.FieldScaling
    for i in range(0, len(rin), 6):
        rin[i+1] /= scaling
        rin[i+3] /= scaling
        rin[i+4] = (rin[i+4] + 1.0 - scaling) / scaling
    return rin
