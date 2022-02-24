# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: __main__.py -- yasim_scripts main frontend.
#
#  VERSION HISTORY:
#  2021-08-15 0.1  : Purposed from the shell script by YU Zhejian.
#  2021-08-21 0.1  : Added; shell script deprecated.
#
# ==============================================================================
"""
__main__.py -- yasim_scripts main frontend.

Please refer to the documentation :doc:`/cmd_interface` to see command-line options.

Description
-----------

This script will:

1) Search for all valid sub-command in :py:mod:`yasim_scripts.main`;
2) Parse arguments provided by user.

If the first non-option (not started by ``-``) is a valid subcommand,
will call the subcommand's :py:func:`main` function with all other arguments and options.
Will raise an error if the subcommand is invalid.

If there's no subcommand,
will display help (if option ``--help`` or ``-h``) or version (with option ``--version`` or ``-v``).

This file also handles log levels.

Subcommand Specifications
-------------------------

A valid subcommand should be a module under :py:mod:`yasim_scripts.main` with :py:func:`main` function defined.
Raw un-parsed arguments EXCEPT sub-module name will be passed to :py:func:`main` function.
"""

import inspect
import logging
import os
import pkgutil
import sys
from typing import List

import yasim_scripts.main
from commonutils.stdlib_helper import logger_helper
from yasim import __version__

__all__ = ['main']

HELP_INFO = """
TODO
"""

lh = logger_helper.get_logger(__name__)

if os.environ.get('LOG_LEVEL') is None:
    logger_helper.set_level('INFO')

valid_subcommand_names = []

for spec in pkgutil.iter_modules(
        yasim_scripts.main.__dict__["__spec__"].submodule_search_locations):
    if not spec.name.startswith("_"):
        valid_subcommand_names.append(spec.name)

_input_subcommand_name = ""
_subcommand_help = "Use 'lscmd' to list all valid subcommands."


def _parse_args(args: List[str]) -> List[str]:
    global _input_subcommand_name, _subcommand_help
    global lh
    HAVE_HELP = False
    HAVE_VERSION = False
    VERBOSE_LEVEL = 0
    i = 0
    while i < len(args):
        name = args[i]
        if name in ('--help', '-h'):
            HAVE_HELP = True
        elif name in ('--__version__', '-v'):
            HAVE_VERSION = True
        elif name in ('--verbose', '-V'):
            VERBOSE_LEVEL += 1
            args.pop(i)
            i -= 1
        elif not name.startswith('-') and _input_subcommand_name == "":
            _input_subcommand_name = name
            args.pop(i)
        i += 1
    if VERBOSE_LEVEL == 1:
        logger_helper.set_level(logging.DEBUG, quiet=False)
    elif VERBOSE_LEVEL >= 2:
        logger_helper.set_level(logger_helper.TRACE, quiet=False)
    if _input_subcommand_name == "lscmd":
        lh.info("Listing modules...")
        for item in valid_subcommand_names:
            print(item)
        sys.exit(0)
    elif _input_subcommand_name == "":
        if HAVE_HELP:
            print(HELP_INFO)
            sys.exit(0)
        elif HAVE_VERSION:
            print(__version__)
            sys.exit(0)
        else:
            lh.exception(f"Subcommand name not set! {_subcommand_help}")
            sys.exit(1)
    return args


def _get_main_func_from_subcommand(name: str):
    """
    Return a subcommands' "main" function.
    """
    global lh
    if name in valid_subcommand_names:
        __import__(f'yasim_scripts.main.{name}')
        i = yasim_scripts.main.__dict__[name]
        if hasattr(i, 'main') and inspect.isfunction(getattr(i, 'main')):
            return i.main


def main():
    """
    Interface accepting raw un-parsed cmdline arguments.

    Should be invoked by :py:mod:`yasim_scripts` only.
    """
    global lh, _input_subcommand_name, _subcommand_help
    if os.environ.get('LOG_LEVEL') is None:
        logger_helper.set_level('INFO')
    lh.info(f'{__doc__.splitlines()[1]} ver. {__version__}')
    lh.info(f'Called by: {" ".join(sys.argv)}')

    args = _parse_args(sys.argv[1:])
    main_fnc = _get_main_func_from_subcommand(_input_subcommand_name)
    if main_fnc is not None:
        sys.exit(main_fnc(args))
    else:
        lh.exception(f"Subcommand '{_input_subcommand_name}' not found! {_subcommand_help}")
        sys.exit(1)


if __name__ == '__main__':
    sys.exit(main())
