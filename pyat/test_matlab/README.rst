Comparing pyAT and AT in Matlab
===============================

It is possible to run both Matlab and Python versions of AT from Python. This
code runs tests to compare the output of the two versions.


Linux Installation & Running the Tests
--------------------------------------

You need an installation of Matlab and a valid licence in order to be able to
call Matlab directly from Python. MATLAB_ROOT is the root directory of your
Matlab installation.

Your Matlab installation will support only certain versions of Python [1]_; find further
information `here <https://uk.mathworks.com/help/matlab/matlab_external/system-
requirements-for-matlab-engine-for-python.html>`_.

Set up a virtualenv using a supported Python version:

* ``cd $AT_ROOT/pyat``
* ``python3 -m venv matlab_venv``
* ``source matlab_venv/bin/activate  # or matlab_venv\Scripts\activate on Windows``
* ``pip install -r requirements.txt``
* ``pip install -e .  # install pyAT into the virtualenv``

Install the Matlab engine for Python, ensuring your virtualenv is still active:

* ``mkdir /tmp/mlp``
* ``cd $MATLAB_ROOT/extern/engines/python``
* ``python setup.py build -b /tmp/mlp install``

For recent Matlab versions, you can also install the Matlab engine for Python
directly from PyPI. You must specify the version of the Matlab engine that matches
your Matlab installation. For example, if you are using Matlab R2025b, you can install
the Matlab engine for Python with:

* ``pip install matlabengine~=25.2.0``

Now run the tests inside your virtualenv:

* ``cd $AT_ROOT/pyat``
* ``python -m pytest test_matlab``


Note: certain versions of GLIBC may be required on Linux: for example,
using R2021a on RHEL7 does not work even though Matlab itself will run.
Using R2021a on RHEL8 does work.


Footnotes
---------

.. [1] `Matlab versions and the Python versions they support: <https://fr.mathworks.com/support/requirements/python-compatibility.html>`_
