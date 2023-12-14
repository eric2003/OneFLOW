# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Scientific Computing'
copyright = '2017~2023, eric'
author = 'eric'
release = '1.0'

show_authors = True

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
#html_theme = 'classic'
#html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

#math_number_all = True  # Set this option to True if you want all displayed math to be numbered. The default is False.
math_eqref_format = 'Eq.{number}'  # gets rendered as, for example, Eq.10.
