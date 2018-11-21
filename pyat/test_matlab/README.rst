Comparing pyAT and AT in Matlab
===============================

It is possible to run both Matlab and Python versions of AT from Python.  This
code runs tests to compare the output of the two versions.


Installation (Linux only)
-------------------------

You need an installation of Matlab and a valid licence in order to be able to
call Matlab directly from Python.  MATLAB_ROOT is the root directory of your
Matlab installation.

Your Matlab installation will support only certain versions of Python, Python3
has proved easiest for me to use.  N.B. please ensure that the versions of the
packages required by pyAT (numpy, scipy and pytest) that you have installed are
compatable with the version of Python that Matlab requires you to run.


Set up a virtualenv:

* ``cd $AT_ROOT``
* ``virtualenv --no-site-packages -p $PYTHON_EXECUTABLE venv``
* ``source venv/bin/activate  # or venv\Scripts\activate on Windows``
* ``pip install -r requirements.txt``
* ``pip install pytest``
* ``python setup.py install``

This should build pyAT into the virtualenv.

Install the Matlab engine for Python into the virtualenv:

* Make sure your virtualenv is still active
* ``mkdir /tmp/mlp``
* ``cd $MATLAB_ROOT/extern/engines/python``
* ``python setup.py build -b /tmp/mlp install``

Now run the tests.

* ``cd $AT_ROOT/pyat``
* ``python -m pytest``
