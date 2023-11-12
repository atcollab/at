# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
# sys.path.insert(0, os.path.abspath('../pyat'))
# print(sys.path)
import at

# -- Project information -----------------------------------------------------

project = 'Accelerator Toolbox'
project_copyright = '2022, ATCollab'
author = 'ATCollab'

# The full version, including alpha/beta/rc tags
release = '.'.join(at.__version__.split('.')[:3])
version = '.'.join(at.__version__.split('.')[:2])


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

root_doc = 'index'
extensions = ['sphinx.ext.autosummary',
              'sphinx.ext.napoleon',
              'sphinx.ext.intersphinx',
              'sphinx.ext.githubpages',
              'sphinx.ext.viewcode',
              'myst_nb',
              'sphinx_copybutton',
              'sphinx_design',
              ]

intersphinx_mapping = {'python': ('https://docs.python.org/3', None),
                       'matplotlib': ('https://matplotlib.org/stable/', None),
                       'numpy': ('https://numpy.org/doc/stable/', None),
                       'scipy': ('https://docs.scipy.org/doc/scipy/', None),
                       }

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["README.rst", "**/*.so", "_build/*"]
rst_prolog = """
.. role:: pycode(code)
   :language: python
"""
autodoc_default_options = {
        # Make sure that any autodoc declarations show the right members
        "members": True,
        "undoc-members": True,
        "inherited-members": False,
        "show-inheritance": True,
        "member-order": "groupwise"
}
autodoc_typehints = 'description'
autodoc_typehints_description_target = 'documented'
autodoc_typehints_format = 'short'
autosummary_generate = True  # Make _autosummary files and include them
autosummary_generate_overwrite = False
autosummary_ignore_module_all = False
autoclass_content = "both"  # include both class docstring and __init__

napoleon_use_rtype = False  # More legible
# napoleon_numpy_docstring = False  # Force consistency, leave only Google
napoleon_custom_sections = [('Returns', 'params_style')]

add_module_names = False

# -- Options for the myst markdown parser ------------------------------------

myst_enable_extensions = [
    "colon_fence",
    "dollarmath",
    "replacements",
    "deflist"
]
myst_heading_anchors = 3
nb_execution_mode = "auto"
nb_execution_allow_errors = True

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
# html_theme = 'sphinx_rtd_theme'
# html_theme = 'sphinx_book_theme'
html_theme = 'pydata_sphinx_theme'
html_logo = 'images/AT.png'
html_copy_source = False
html_theme_options = {
    "github_url": "https://github.com/atcollab/at",
}
html_sidebars = {
    "index": [],
    "common/about": [],
}
# creates an additional page, but impossible to link to it...
# if os.environ.get('READTHEDOCS') == 'True':
#     html_additional_pages = {
#         "index": "plink.html"
# }

html_css_files = ["css/custom_at.css"]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static', 'atdocs']


# -- Options for copybutton  -------------------------------------------------

copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True
copybutton_only_copy_prompt_lines = True
