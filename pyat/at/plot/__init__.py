"""AT plotting functions"""
print('plot start')
try:
    import matplotlib
    from .synopt import *
    from .generic import *
    from .specific import *
except (ImportError, RuntimeError) as exc:
    print(exc)
    print('Plotting is disabled')
print('plot end')
