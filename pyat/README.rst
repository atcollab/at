pyAT
====

pyAT is a Python interface to the pass methods defined in Accelerator Toolbox,
implemented by compiling the C code used in the AT 'integrators' plus a Python
extension.

It supports Python 2.7 and 3.3 to 3.6.


Installation preparation (Windows)
----------------------------------

Download Microsoft Visual C++ Compiler for Python 2.7 (`here
<https://www.microsoft.com/en-us/download/details.aspx?id=44266>`_), and use
the Visual C++ Command Prompt of the correct architecture to build pyat.

For newer versions of Python you need the appropriate version of Visual C++.


Installation (all platforms)
----------------------------

All the binaries should be built when building the Python extension.

It is easiest to do this using a virtualenv:

* ``virtualenv --no-site-packages venv``
* ``source venv/bin/activate  # or venv\Scripts\activate on Windows``
* ``pip install numpy``
* ``pip install scipy``
* ``pip install pytest``
* ``python setup.py develop  # install into the virtualenv``

Finally, you should be able to run the tests:

* ``py.test test``

Any changes to .py files are automatically reinstalled in the build, but to
ensure any changes to .c files are reinstalled rerun the following:

* ``python setup.py install``


Comparing results with Matlab
-----------------------------

There is a second set of tests that require a Matlab licence and allows
comparing results directly with a Matlab session.  See test_matlab/README
for information.


Debugging
---------

Print statements in the C code will work once the integrators are
recompiled.  To force recompilation, remove the build directory:

* ``rm -rf build``
