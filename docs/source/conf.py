# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import datetime
import os
import sys

sys.path.insert(0, os.path.abspath('../../'))
import qcdb

# -- Project information -----------------------------------------------------

project = 'QCDB'
copyright = f'2013-{datetime.datetime.today().year}, The {project} Project'
author = 'The QCDB Development Team'

# The short X.Y version
version = qcdb.__version__
# The full version, including alpha/beta/rc tags
release = qcdb.__version__

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
needs_sphinx = '3.5'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # from Sphinx
    'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.extlinks',
    'sphinx.ext.graphviz',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.graphviz',
    "sphinx_autodoc_typehints",
    "sphinx.ext.githubpages",
    # from Astropy
    'sphinx_automodapi.automodapi',
    'sphinx_automodapi.automodsumm',
    'sphinx_automodapi.smart_resolver',
    # from Cloud
    'cloud_sptheme.ext.index_styling',
    'cloud_sptheme.ext.escaped_samp_literals',
    # from Psi4
    'sphinx_psi_theme.ext.psidomain',
    'sphinx_psi_theme.ext.relbar_toc',
]

autosummary_generate = True
automodapi_toctreedirnm = 'api'
#numpydoc_show_class_members = False
#automodsumm_inherited_members = True
autodoc_typehints = "description"
napoleon_use_param = True
napoleon_use_rtype = True





# Import Sphinx themes
import sphinx_psi_theme


autodoc_default_flags = ['members',
                         'undoc-members',
                         'inherited-members',  # disabled because there's a bug in sphinx
                         'show-inheritance',
                        ]




# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# Suppress warnings (sphinx 1.4.2)
suppress_warnings = ['image.nonlocal_uri']

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['abbr_accents.rst']

# The reST default role (used for this markup: `text`) to use for all documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
show_authors = False
#show_authors = True

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'default'

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []


# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_psi_theme'
import sphinx_psi_theme

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {"roottarget": "index",
                      "barlink_github": "http://github.com/qcdb/qcdb",
                      "barlink_home": "https://github.com/qcdb/qcdb/blob/master/README.md",
                      "barlink_help": "http://forum.psicode.org",
                      "barlink_docsrc": "docs/sphinx",
                      "barlink_pkgnamehtml": """<span style="font-family: Optima, sans-serif;">QCDB</span>""",
                      "footerbgcolor": '#555465',
                      "relbarbgcolor": '#6a7591',
                     }

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = [sphinx_psi_theme.get_theme_dir()]

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = ''

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "_static/qcdbsquare.png"

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = '_static/favicon-qcdb.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%A, %d %B %Y %I:%M%p'

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {'**': ['localtoc.html', 'relations.html', 'sourcelink.html', 'searchbox.html']}
html_sidebars = {'**': ['searchbox.html', 'localtoc.html'],
                 'index': ['searchbox.html', 'impressum.html']}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
#html_show_sourcelink = True

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'qcdbdoc'

intersphinx_mapping = {'python': ('http://docs.python.org/3.6', None),
                       'numpy': ('http://docs.scipy.org/doc/numpy/', None),
                       'scipy': ('http://docs.scipy.org/doc/scipy/reference/', None),
                       'matplotlib': ('http://matplotlib.sourceforge.net/', None),
                       'psi4.core': ('http://psicode.org/psi4manual/master/', None),
                      }

nbsphinx_allow_errors = True
nbsphinx_execute = 'never'
nbsphinx_timeout = 180

# Abbreviations
rst_epilog = """

.. color: #273896;">
.. |PSIfour| raw:: html

    <span style="font-family: Optima, sans-serif; text-transform: none;">P<span style="font-size: 82%;">SI</span>4</span>

.. |PSIfours| replace:: |PSIfour|\ 's
.. |qcdb| raw:: html

    <span style="font-family: Optima, sans-serif; text-transform: none;">QCDB</span>

.. |qcdbs| replace:: |qcdb|\ 's
.. |dl| replace:: :math:`\Rightarrow`
.. |dr| replace:: :math:`\Leftarrow`
.. |kcalpermol| replace:: kcal mol\ :sup:`-1`
.. |Angstrom| replace:: |AA|\ ngstr\ |o_dots|\ m
.. |MollerPlesset| replace:: M\ |o_slash|\ ller--Plesset
.. |--| unicode:: U+02013 .. en dash
   :trim:
.. |w--w| unicode:: U+02013 .. en dash
.. |---| unicode:: U+02014 .. em dash
   :trim:
.. |w---w| unicode:: U+02014 .. em dash
.. include:: /abbr_accents.rst
"""

rst_prolog = """
.. highlight:: python
   :linenothreshold: 1
"""

extlinks = {'source':    ('https://github.com/qcdb/qcdb/blob/master/%s', 'qcdb/'),
            #'srcsample': ('https://github.com/qcdb/qcdb/blob/master/samples/%s/input.dat', ''),
            'srcbasis':  ('https://github.com/qcdb/qcdb/blob/master/share/qcdb/basis/%s.gbs', ''),
            #'srcplugin': ('https://github.com/qcdb/qcdb/blob/master/plugins/%s', ''),
            #'srcefpfrag':('https://github.com/ilyak/libefp/blob/master/fraglib/%s.efp', ''),
            'srcdb':     ('https://github.com/qcdb/qcdb/blob/master/share/qcdb/databases/%s.py', '') }

def setup(app):

    app.add_object_type('qcvar', 'qcvar', indextemplate='single: %s')
