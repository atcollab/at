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

Now run the tests inside your virtualenv:

* ``cd $AT_ROOT/pyat``
* ``python -m pytest test_matlab``


Note: certain versions of GLIBC may be required on Linux: for example,
using R2021a on RHEL7 does not work even though Matlab itself will run.
Using R2021a on RHEL8 does work.


Footnotes
---------

.. [1] Matlab versions and the Python versions they support:

   +----------------+--------------------------+
   | Matlab Release | Supported Python Version |
   +================+==========================+
   |  2021b         |  2.7, 3.7, 3.8, 3.9      |
   +----------------+--------------------------+
   |  2021a         |  2.7, 3.7, 3.8           |
   +----------------+--------------------------+
   |  2020b         |  2.7, 3.6, 3.7, 3.8      |
   +----------------+--------------------------+
   |  2020a         |  2.7, 3.6, 3.7           |
   +----------------+--------------------------+
   |  2019b         |  2.7, 3.6, 3.7           |
   +----------------+--------------------------+
   |  2019a         |  2.7, 3.5, 3.6, 3.7      |
   +----------------+--------------------------+
   |  2018b         |  2.7, 3.5, 3.6           |
   +----------------+--------------------------+
   |  2018a         |  2.7, 3.5, 3.6           |
   +----------------+--------------------------+
   |  2017b         |  2.7, 3.4, 3.5, 3.6      |
   +----------------+--------------------------+
   |  2017a         |  2.7, 3.4, 3.5           |
   +----------------+--------------------------+
   |  2016b         |  2.7, 3.3, 3.4, 3.5      |
   +----------------+--------------------------+
   |  2016a         |  2.7, 3.3, 3.4           |
   +----------------+--------------------------+
   |  2015b         |  2.7, 3.3, 3.4           |
   +----------------+--------------------------+
   |  2015a         |  2.7, 3.3, 3.4           |
   +----------------+--------------------------+
   |  2014b         |  2.7, 3.3                |
   +----------------+--------------------------+
   |  <=2014a       |  Not supported           |
   +----------------+--------------------------+

