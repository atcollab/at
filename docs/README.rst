##################
PyAT Documentation
##################

This directory is the source of the PyAT Documentation

**********************************
Contributing to PyAT documentation
**********************************

Documentation tools
===================

The PyAT documentation is compiled with
`Sphinx <https://www.sphinx-doc.org/en/master/index.html>`_.

The documentation on GitHub Pages is compiled at each ``push`` action on the
``master`` branch of AT. This ensures that the documentation is immediately
updated for any modification in python docstrings.

For a local check when writing documentation, Sphinx and the
necessary extensions are automatically installed if you install PyAT with
the [doc] option::

    $ pip install -e ".[doc]"

Then, the compilation may be triggered from the ``docs`` directory with::

    $ cd <at_root>/docs
    $ make html

Writing documentation
=====================

The API documentation is extracted from the python docstrings. PyAT uses the
`Google style guide <https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings>`_
for docstrings. Sphinx is rather picky with indentation and blank lines, so
better looking at a working example while writing.

Any other documentation may be added by using
`reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_
or `markdown <https://myst-parser.readthedocs.io/en/latest/syntax/syntax.html>`_,
and putting files in the ``<at_root>/docs`` directory or any subdirectory.
Links to new files must be provided in ``<at_root>/docsrc/index.rst``, or in
cross-references from another file.