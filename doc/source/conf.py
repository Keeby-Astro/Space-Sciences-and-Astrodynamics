# Sphinx configuration for Space-Sciences-and-Astrodynamics
import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../curtis_scripts/python')))

project = 'Space-Sciences-and-Astrodynamics'
author = 'Keeby Astro'
release = '0.1.0a0'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx_autodoc_typehints',
]
html_theme = 'sphinx_rtd_theme'

# -- General information ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

# The master toctree document.
master_doc = 'index'

# -- Options for autodoc ---------------------------------------------------
autodoc_member_order = 'bysource'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}

# -- Options for napoleon --------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# -- Paths for static and templates ----------------------------------------
html_static_path = []
templates_path = []


def preserve_script_docstrings(app, what, name, obj, options, lines):
    # The educational parameter descriptions use aligned, fixed-width text.
    if lines:
        lines[:] = ['::', ''] + ['    ' + line for line in lines]


def setup(app):
    app.connect('autodoc-process-docstring', preserve_script_docstrings, priority=100)