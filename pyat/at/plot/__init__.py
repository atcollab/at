"""AT plotting functions"""
try:
    # noinspection PyPackageRequirements
    import matplotlib
    from .synopt import *
    from .generic import *
    from .specific import *
    from .standalone import *
except (ImportError, RuntimeError) as exc:
    print(exc)
    print('Plotting is disabled')
