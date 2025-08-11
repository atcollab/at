"""AT plotting functions"""

try:
    # noinspection PyPackageRequirements
    import matplotlib
except (ImportError, RuntimeError) as exc:
    print("matplotlib is unavailable => plotting is disabled.")
    print("To enable plotting functions, run \"pip install matplotlib\".")
else:
    from .synopt import *
    from .generic import *
    from .specific import *
    from .standalone import *
    from .resonances import *
    from .response_matrix import *
