pyAT Developer Notes
====================

See also README.rst.


Installation preparation (Windows)
----------------------------------

Download Microsoft Visual C++ Compiler for Python 2.7 (`here
<https://www.microsoft.com/en-us/download/details.aspx?id=44266>`_), and use
the Visual C++ Command Prompt of the correct architecture to build pyat.

For newer versions of Python you need the appropriate version of Visual C++.


Installation (all platforms)
----------------------------

All the binaries should be built when building the Python extension.

It is easiest to do this using a virtualenv. inside pyat:

Python 2 support will be removed soon. If you are still using Python 2, you
need virtualenv installed:

* ``virtualenv --no-site-packages venv``

If you are using Python 3, you can use the built-in venv module:

* ``python3 -m venv venv``

Then:

* ``source venv/bin/activate  # or venv\Scripts\activate on Windows``
* ``pip install -r requirements.txt``
* ``pip install -e .``

Finally, you should be able to run the tests:

* ``python -m pytest test``


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

Any changes to .py files are automatically reinstalled in the build, but to
ensure any changes to .c files are reinstalled rerun:

* ``python setup.py develop``

If you get strange behaviour even after running setup.py develop again, then
running the following, inside pyat, should fix it:

* ``rm -rf build``
* ``find at -name "*.pyc" -exec rm '{}' \;``
* ``find at -name "*.so" -exec rm '{}' \;``
* ``python setup.py develop``

N.B. setup.py develop needs to be run with the same version of Python (and
numpy) that you are using to run pyAT.

Releasing a version to PyPI
---------------------------

Because pyAT compiles C code, releasing a version is not simple. The code
must be compiled for different operating systems and Python versions.

To do this, we use the continuous integration service Github Actions.
When a tag of the form pyat-x.y.z is pushed to Github, wheels for each
supported platform will be built and automatically uploaded as an 'artifact'.

Once there, someone should manually test that the wheels are working correctly,
then they can manually download the files and upload them to PyPI itself.
Note that this was 46 different files for pyat-0.0.4 covering different platforms and
architectures.

The configuration for this is in .github/workflows/build-wheels.yml.
