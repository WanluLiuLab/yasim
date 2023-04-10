"""
Configuration file for the Sphinx documentation builder.
"""

# pylint: disable=wrong-import-position, invalid-name

import glob
import os
import shutil
from collections import defaultdict

import tomli

import yasim
from labw_utils.commonutils import libfrontend
from labw_utils.devutils.sphinx_helper import convert_ipynb_to_myst

os.environ['SPHINX_BUILD'] = '1'  # Disable chronolog.
THIS_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(THIS_DIR)


def generate_cli_docs(
        config_toml_file_path: str,
        dest_dir_path: str
):
    os.makedirs(dest_dir_path, exist_ok=True)
    with open(config_toml_file_path, "rb") as toml_reader:
        config_toml = tomli.load(toml_reader)

    arg_parsers = defaultdict(lambda: [])
    for main_module in config_toml["names"]:
        for subcommand in libfrontend.get_subcommands(main_module):
            parser = libfrontend.get_argparser_from_subcommand(main_module, subcommand)
            this_help_path = os.path.join(dest_dir_path, f"{main_module}.{subcommand}.txt")
            if parser is not None:
                with open(this_help_path, "w") as writer:
                    writer.write(parser.format_help())
                arg_parsers[main_module].append(subcommand)
            else:
                doc = libfrontend.get_doc_from_subcommand(main_module, subcommand)
                if doc is None:
                    continue
                else:
                    # doc_sio = io.StringIO(doc)
                    with open(this_help_path, "w") as writer:
                        writer.write(doc)
                    arg_parsers[main_module].append(subcommand)

    with open(os.path.join(dest_dir_path, "index.md"), "w") as index_writer:
        index_writer.write("# Command-Line Interfaces\n\n")
        for main_module, subcommands in arg_parsers.items():
            main_module_correct_name = main_module.replace("._main", "").replace(".main", "")
            index_writer.write(f"## `{main_module_correct_name}`\n\n")
            for subcommand in subcommands:
                index_writer.write(f"### `{subcommand}`\n\n")
                index_writer.write(
                    "```{literalinclude} " + f"{main_module}.{subcommand}.txt" + "\n:language: text\n```\n\n"
                )


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
convert_ipynb_to_myst(THIS_DIR)

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
    "*.ipynb"
]
myst_enable_extensions = ["deflist", "dollarmath"]
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
    ".ipynb.md": ["jupytext.reads", {"fmt": "md:myst"}]
}

# Insert both docstring of the class and constructor.
autodoc_default_options = {
    'special-members': '__init__',
}

# Intersphinx settings
intersphinx_mapping = {
    'python': ('https://docs.python.org/3.8', None),
    'joblib': ('https://joblib.readthedocs.io/en/latest', None),
}

# myst-nb settings
nb_execution_timeout = 1200
nb_execution_mode = "cache"
nb_merge_streams = True
nb_execution_allow_errors = True

generate_cli_docs(
    os.path.join(THIS_DIR, "cli_docs.toml"),
    os.path.join(THIS_DIR, "_cli_docs")
)
