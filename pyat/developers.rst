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

We use `Semantic Versioning <https://semver.org/>`_ for this project. As documented
by the full specification, we use a version number with three dot-separated
numbers: ``MAJOR.MINOR.PATCH``. Increment the:

* MAJOR version when you make incompatible API changes
* MINOR version when you add functionality in a backward compatible manner
* PATCH version when you make backward compatible bug fixes

Release procedure
-----------------

To upload a release to PyPI, you will need to be a 'maintainer' of
`Accelerator Toolbox on PyPI <https://pypi.org/project/accelerator-toolbox/>`_.

For testing any version that you have installed, the simple snippet in
``README.rst`` is sufficient.

Decide the Python versions that should be supported in the release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Set these Python versions in ``requires-python`` in ``pyproject.toml``
* Set at least these Python versions as ``python-version`` in ``.github/workflows/python-tests.yml``

Determine the minimum ``numpy`` and ``scipy`` versions:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* the version ensuring the requirements necessary to **run** PyAT is set in the
  ``dependencies`` item of the ``[project]`` section of ``pyproject.toml``
* The version required to **build** PyAT is set in the ``requires`` item of the
  ``[build-system]`` section of ``pyproject.toml``. It depends on the python version
  and must be higher or equal to the "run" version.
* To avoid ABI compatibility issues, the pre-compiled binaries are built with the
  earliest possible version of numpy for the given Python version. This ensures that
  the user's libraries are more recent than the one AT has been compiled with. For
  that, a copy of ``pyproject.toml`` named ``githubproject.toml`` is used for
  compilation. In this copy, the numpy version specifications are set using ``~=``
  instead of minimum (``>=``). Apart from these lines, the 2 files
  should be strictly identical.

Prepare the "Release notes"
~~~~~~~~~~~~~~~~~~~~~~~~~~~
A draft can be obtained by creating a new tag on GitHub. Click "Draft a new release"
in the releases page, choose a new tag in the form ``pyat-x.y.z`` with the correct
incremented version number. The ``pyat-`` prefix is necessary to identify python releases.
Select the master branch as target. In the description area, choose the current
release in the "previous tag" pull-down, and press "Generate release notes".

The generated notes can now be copied and edited. You can then either cancel or
save the release as a draft while editing the release notes.

The ``## What's changed`` section should be split into ``## Bug fixes`` and
``## New features``. It must be filtered to keep only the python changes, ignoring
the Matlab ones. The tags on each pull request are there to help in this filtering.

The release notes should start with a ``## Main modifications`` section summarising
the important points of the new release.

They must end with a section pointing out ``## Incompatibilities`` and mentioning the
necessary actions before upgrading to this release.

Open a Pull Request for the new release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The goal is to make all contributors aware of the new release, to check that no pending
modifications are worth being included and to review the release notes.

There should be no code modifications except updates of version informations in the
documentation. Once the pull request is approved and merged, the release may be built.


Build the release
~~~~~~~~~~~~~~~~~

Now either go back to the draft release saved above, or start again the procedure,
but now finalising with the ``Publish`` button.

If all goes well, there will be a build of "Build and upload wheels and sdist"
associated with the tag ``pyat-x.y.z``: on the `Github Actions page <https://github.com/atcollab/at/actions/workflows/build-python-wheels.yml>`_. This build will have
'tar.gz' and 'wheels' downloads available.

Upload the release to ``pip``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Download the tar.gz and wheels files and unzip them into a directory ``<dir>``
* Manually install at least one wheel to make sure that it has built correctly
* Install Twine for uploading the files to PyPI. One way is to use
  `pipx <https://pypa.github.io/pipx/>`_. Once pipx installed, use::

     $ pipx install twine  # or:
     $ pipx upgrade twine

* Use Twine to upload the files to PyPI. You will be prompted for your PyPI credentials::

    $ twine upload <dir>/*.whl
    $ twine upload <dir>/*.tar.gz

* Finally, check that the wheels are uploaded properly. You can use the same virtualenv::

    $ pip install accelerator-toolbox

Note that 46 different files were uploaded for pyat-0.0.4 covering different
platforms and architectures.

The configuration for this is in ``.github/workflows/build-python-wheels.yml``.
