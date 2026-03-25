"""Sphinx configuration for Threads documentation."""

import os
import sys

# Add project src to path so autodoc can import threads_arg
sys.path.insert(0, os.path.abspath(os.path.join('..', 'src')))

project = 'Threads'
author = 'Threads authors'
copyright = '2024-2025, Threads authors'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx_autodoc_typehints',
]

# Napoleon settings (Google/NumPy style docstrings)
napoleon_google_docstrings = True
napoleon_numpy_docstrings = True
napoleon_include_init_with_doc = True

# Autodoc settings
autodoc_member_order = 'bysource'
autodoc_typehints = 'description'
autosummary_generate = True

# Intersphinx mappings
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}

# Theme
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'navigation_depth': 3,
    'collapse_navigation': False,
}

# General settings
exclude_patterns = ['_build']
templates_path = ['_templates']
html_static_path = ['_static']

os.makedirs(os.path.join(os.path.dirname(__file__), '_static'), exist_ok=True)
