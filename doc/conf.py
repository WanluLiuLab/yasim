"""
Configuration file for the Sphinx documentation builder.
"""

# pylint: disable=wrong-import-position, invalid-name

import glob
import os
import shutil
import sys

import tomli

os.environ['SPHINX_BUILD'] = '1'  # Disable chronolog.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(THIS_DIR)

# Enable scan of packages in src which is not installed.
sys.path.insert(0, os.path.join(ROOT_DIR, 'src'))

import yasim


def scan_dir(path_to_scan: str):
    """
    Recursively list files in directories.

    <https://www.sethserver.com/python/recursively-list-files.html>
    """

    files = []
    dirlist = [path_to_scan]
    while len(dirlist) > 0:
        for (dirpath, dirnames, filenames) in os.walk(dirlist.pop()):
            dirlist.extend(dirnames)
            files.extend(map(lambda n: os.path.join(*n), zip([dirpath] * len(filenames), filenames)))
    return files


def copy_doc_files(from_path: str, to_path: str):
    """
    Copy items to project root
    """
    os.makedirs(to_path, exist_ok=True)
    for name in glob.glob(from_path):
        shutil.copy(name, to_path + os.sep)


copy_doc_files(os.path.join(ROOT_DIR, '*.md'), os.path.join(THIS_DIR, "_root"))

# -- Project information -----------------------------------------------------

with open(os.path.join(ROOT_DIR, "pyproject.toml"), "rb") as reader:
    parsed_pyproject = tomli.load(reader)

project = parsed_pyproject["project"]["name"]
author = "&".join([author["name"] for author in parsed_pyproject["project"]["authors"]])
copyright_string = f'2022-2023, {author}'
release = yasim.__version__

# -- General configuration ---------------------------------------------------

html_theme = 'sphinx_rtd_theme'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.todo',
    #    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    "sphinx.ext.viewcode",
    'myst_nb',
    'sphinx_copybutton',
    'sphinxcontrib.bibtex'
]
myst_enable_extensions = ["deflist"]
bibtex_bibfiles = ['refs.bib']
exclude_patterns = [
    '_build',
    'Thumbs.db',
    '.DS_Store',
    '.virtualenv/**'
]

# html_static_path = ['_static']

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'myst-nb',
}

nb_custom_formats = {
    ".ipynb.py": ["jupytext.reads", {"fmt": "py:percent"}]
}

# Insert both docstring of the class and constructor.
autodoc_default_options = {
    'special-members': '__init__',
}

# Intersphinx settings
intersphinx_mapping = {
    'python': ('https://docs.python.org/3.8', None),
    'joblib': ('https://joblib.readthedocs.io/en/latest', None),
    'sklearn': ('https://scikit-learn.org/stable', None),
    'torch': ('https://pytorch.org/docs/stable', None),
}

# myst-nb settings
nb_execution_timeout = 1200
nb_execution_mode = "cache"
nb_merge_streams = True
