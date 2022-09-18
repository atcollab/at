pyAT Developer Notes
====================

See also README.rst.


Installation preparation
------------------------

Windows
~~~~~~~

Download the appropriate version of Visual C++.

Linux
~~~~~

Install the Python development package for your OS. For example, using yum::

    $ sudo yum install python3-devel

Source download
---------------
Download the latest version of AT::

    $ git clone https://github.com/atcollab/at.git

Installation (all platforms)
----------------------------

All the binaries should be built when building the Python extension.

It is easiest to do this using a virtualenv. inside pyat::

    $ python3 -m venv venv

Then:

* activate the virtual environment::

    $ source venv/bin/activate  # or venv\Scripts\activate on Windows

* make sure you have a recent pip installer::

    $ pip install --upgrade pip

* Go to the AT root directory::

    cd <atroot>

* install AT, with the developer tools.
  2 option sets are available:

  * ``[dev]`` installs test modules (``pytest``,...)
  * ``[doc]`` installs documentation tools (``sphinx`` and extensions)

  They can be installed together: [dev, doc]::

    $ pip install -e ".[dev]"

Finally, you should be able to run the tests::

    $ python -m pytest pyat/test


Comparing results with Matlab
-----------------------------

There is a second set of tests that require a Matlab licence and allows
comparing results directly with a Matlab session.  See test_matlab/README
for information.


Debugging
---------

Print statements in the C code will work once the integrators are
recompiled.  To force recompilation, remove the build directory::

    $ rm -rf build

Any changes to .py files are automatically reinstalled in the build, but to
ensure any changes to .c files are reinstalled rerun::

    $ pip install -e .

If you get strange behaviour even after running pip install develop again, then
running the following, inside pyat, should fix it::

    $ rm -rf build
    $ find at -name "*.pyc" -exec rm '{}' \;
    $ find at -name "*.so" -exec rm '{}' \;
    $ pip install -e .

N.B. ``pip install -e .`` needs to be run with the same version of Python (and
numpy) that you are using to run pyAT.

Releasing a version to PyPI
---------------------------

Because pyAT compiles C code, releasing a version is not simple. The code
must be compiled for different operating systems and Python versions.

To do this, we use the continuous integration service Github Actions.
When a tag of the form pyat-* is pushed to Github, wheels for each
supported platform will be built and automatically uploaded as an 'artifact'.

Release procedure
-----------------

To upload a release to PyPI, you will need to be a 'maintainer' of
`Accelerator Toolbox on PyPI <https://pypi.org/project/accelerator-toolbox/>`_.

For testing any version that you have installed, the simple snippet in
``README.rst`` is sufficient.

* Decide the Python versions that should be supported in the release
   * Set these Python versions in ``python_requires`` in ``setup.cfg``
   * Set at least these Python versions as ``python-version`` in ``.github/workflows/python-tests.yml``
* Determine the minimum Numpy version that is required for those Python versions
   * Set this numpy version in ``install_requires`` in ``setup.cfg``
* Push a tag ``pyat-x.y.z`` to Github

If all goes well, there will be a build of "Build and upload wheels and sdist"
associated with the tag ``pyat-x.y.z``: on the `Github Actions page <https://github.com/atcollab/at/actions/workflows/build-python-wheels.yml>`_. This build will have
'tar.gz' and 'wheels' downloads available.

* Download the tar.gz and wheels files and unzip them into a directory ``<dir>``
* Manually install at least one wheel to make sure that it has built correctly
* Install Twine for uploading the files to PyPI. One way to do this is to create a new virtualenv::

    $ python3 -m venv venv
    $ source venv/bin/activate
    $ pip install twine

* Use Twine to upload the files to PyPI. You will be prompted for your PyPI credentials::

    $ twine upload <dir>/*.whl
    $ twine upload <dir>/*.tar.gz

* Finally, check that the wheels are uploaded properly. You can use the same virtualenv::

    $ pip install accelerator-toolbox

Note that 46 different files were uploaded for pyat-0.0.4 covering different
platforms and architectures.

The configuration for this is in .github/workflows/build-python-wheels.yml.
