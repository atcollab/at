import os
import sys
try:
    import matlab.engine
except ImportError:
    print('Matlab comparison tests require Matlab Python Engine installed.')
    print('Python will exit.')
    sys.exit()


ROOT_DIR = os.path.realpath(os.path.join(os.path.dirname(__file__), '../..'))
LATTICE = os.path.join(ROOT_DIR, 'atmat/atdemos/atmatchExamples/ExampleATMATCH/dba.mat')


def initialise_matlab():
    eng = matlab.engine.start_matlab()
    eng.addpath(eng.genpath(os.path.join(ROOT_DIR, 'atintegrators/')))
    eng.addpath(eng.genpath(os.path.join(ROOT_DIR, 'atmat/')))
    return eng
