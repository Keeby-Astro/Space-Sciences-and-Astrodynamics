# Sphinx configuration for Space-Sciences-and-Astrodynamics
import os
import sys
sys.path.insert(0, os.path.abspath('../../curtis_scripts/python'))

project = 'Space-Sciences-and-Astrodynamics'
author = 'Keeby Astro'
release = '0.1.0a0'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc_typehints',
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
html_static_path = ['_static']
templates_path = ['_templates']
