Comparing pyAT and AT in Matlab
===============================

It is possible to run both Matlab and Python versions of AT from Python. This
code runs tests to compare the output of the two versions.


Linux Installation & Running the Tests
--------------------------------------

You need an installation of Matlab and a valid licence in order to be able to
call Matlab directly from Python. MATLAB_ROOT is the root directory of your
Matlab installation.

Your Matlab installation will support only certain versions of Python [1]_; if
you're using Python 2 it needs to be compiled using ucs4 [2]_, find further
information `here <https://uk.mathworks.com/help/matlab/matlab_external/system-
requirements-for-matlab-engine-for-python.html>`_. N.B. please ensure that the
versions of the packages required by pyAT (numpy, scipy and pytest) that you
have installed are compatible with the version of Python that Matlab requires
you to run.

Set up a virtualenv using a supported Python version:

* ``cd $AT_ROOT``
* ``virtualenv --no-site-packages -p $PYTHON_EXECUTABLE venv``
* ``source venv/bin/activate  # or venv\Scripts\activate on Windows``
* ``pip install -r requirements.txt``
* ``python setup.py install  # install pyAT into the virtualenv``

Install the Matlab engine for Python, ensuring your virtualenv is still active:

* ``mkdir /tmp/mlp``
* ``cd $MATLAB_ROOT/extern/engines/python``
* ``python setup.py build -b /tmp/mlp install``

Now run the tests inside your virtualenv:

* ``cd $AT_ROOT/pyat``
* ``$PYTHON_EXECUTABLE -m pytest test_matlab``


Footnotes
---------

.. [1] Matlab versions and the Python versions they support:

   +----------------+--------------------------+
   | Matlab Release | Supported Python Version |
   +================+==========================+
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

.. [2] To check if your Python version is compiled with ucs2 or ucs4::

   >>> import sys
   >>> print('ucs4' if sys.maxunicode > 65535 else 'ucs2')

