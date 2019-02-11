import os
import sys
try:
    import matlab.engine
except ImportError:
    print('Matlab comparison tests require Matlab Python Engine installed.')
    print('Python will exit.')
    sys.exit()


ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '../..'))
dba_ring = os.path.join(ROOT_DIR, 'pyat/test_matlab/dba.mat')
hmba_ring = os.path.join(ROOT_DIR, 'pyat/test_matlab/hmba.mat')


def initialise_matlab():
    eng = matlab.engine.start_matlab()
    eng.addpath(eng.genpath(os.path.join(ROOT_DIR, 'atintegrators/')))
    eng.addpath(eng.genpath(os.path.join(ROOT_DIR, 'atmat/')))
    eng.addpath(os.path.dirname(__file__))
    return eng
