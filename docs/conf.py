# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "argweavers"
copyright = "2013-2023, Cornell University; 2024, Tang Ziya"
author = "Tang Ziya"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "breathe",
    "myst_parser",
    "autodoc2",
    "sphinx.ext.graphviz",
    "sphinx.ext.napoleon",
]

autodoc2_packages = [
    "../argweaver",
]
autodoc2_render_plugin = "myst"
breathe_projects = {"argweaver": "../"}
graphviz_output_format = "svg"

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"
html_static_path = ["_static"]
