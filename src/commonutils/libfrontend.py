import importlib
import inspect
import logging
import os
import pkgutil
import sys
from typing import List, Iterable, Callable

from commonutils.stdlib_helper import logger_helper

__all__ = ['setup_frontend']

lh = logger_helper.get_logger(__name__)

if os.environ.get('LOG_LEVEL') is None:
    logger_helper.set_level('INFO')


def _get_subcommands(package_main_name: str) -> Iterable[str]:
    for spec in pkgutil.iter_modules(pkgutil.resolve_name(package_main_name).__spec__.submodule_search_locations):
        if not spec.name.startswith("_"):
            yield spec.name


def _get_main_func_from_subcommand(
        package_main_name: str,
        subcommand_name: str
) -> Callable[[List[str]], int]:
    """
    Return a subcommands' "main" function.
    """
    importlib.import_module(f'{package_main_name}.{subcommand_name}')
    i = pkgutil.resolve_name(f'{package_main_name}.{subcommand_name}')
    if hasattr(i, 'main') and inspect.isfunction(getattr(i, 'main')):
        return i.main


def lscmd(
        valid_subcommand_names: Iterable[str]
):
    lh.info("Listing modules...")
    for item in valid_subcommand_names:
        print(item)
    sys.exit(0)


class _ParsedArgs:
    input_subcommand_name: str = ""
    have_help: bool = False
    have_version: bool = False
    verbose_level: int = 0
    parsed_args: List[str] = []

    def set_verbose_level(self):
        if os.environ.get('LOG_LEVEL') is not None:
            return
        if self.verbose_level == 1:
            logger_helper.set_level(logging.DEBUG, quiet=False)
        elif self.verbose_level >= 2:
            logger_helper.set_level(logger_helper.TRACE, quiet=False)


def _parse_args(
        args: List[str]
) -> _ParsedArgs:
    parsed_args = _ParsedArgs()
    i = 0
    while i < len(args):
        name = args[i]
        if name in ('--help', '-h'):
            parsed_args.have_help = True
        elif name in ('--version', '-v'):
            parsed_args.have_version = True
        elif name in ('--verbose', '-V'):
            parsed_args.verbose_level += 1
            args.pop(i)
            i -= 1
        elif not name.startswith('-') and parsed_args.input_subcommand_name == "":
            parsed_args.input_subcommand_name = name
            args.pop(i)
        i += 1
    parsed_args.parsed_args = args
    return parsed_args


def _format_help_info(package_main_name:str) -> str:
    return f"""
This is frontend of `{package_main_name.split('.')[0].strip()}` provided by `commonutils.libfrontend`.

SYNOPSYS: {sys.argv[0]} [[SUBCOMMAND] [ARGS_OF SUBCOMMAND] ...] [-h|--help] [-v|--version] [-V|--verbose]

If a valid [SUBCOMMAND] is present, will execute [SUBCOMMAND] with all other arguments
    except [-V|--verbose], which is used to increase log levels.

If no valid [SUBCOMMAND] is present, will fail to errors.

If no [SUBCOMMAND] is present, will consider options like:
    [-h|--help] show this help message and exit
    [-v|--version] show package version and other information

Use `lscmd` as subcommand with no options to see available subcommands.
"""

def _act_on_args(
        parsed_args: _ParsedArgs,
        package_main_name: str,
        version: str,
        help_info: str,
        subcommand_help: str

):
    parsed_args.set_verbose_level()
    valid_subcommand_names = _get_subcommands(package_main_name)
    if parsed_args.input_subcommand_name == "lscmd":
        lscmd(valid_subcommand_names)
    elif parsed_args.input_subcommand_name == "":
        if parsed_args.have_help:
            if help_info is None:
                help_info = _format_help_info(package_main_name)
            print(help_info)
            sys.exit(0)
        elif parsed_args.have_version:
            print(version)
            sys.exit(0)
        else:
            lh.exception(f"Subcommand name not set! {subcommand_help}")
            sys.exit(1)
    elif parsed_args.input_subcommand_name in valid_subcommand_names:
        main_fnc = _get_main_func_from_subcommand(
            package_main_name=package_main_name,
            subcommand_name=parsed_args.input_subcommand_name
        )
        if main_fnc is not None:
            sys.exit(main_fnc(parsed_args.parsed_args))
        else:
            lh.exception(f"Subcommand '{parsed_args.input_subcommand_name}' not found! {subcommand_help}")
            sys.exit(1)
    else:
        lh.exception(f"Subcommand '{parsed_args.input_subcommand_name}' not found! {subcommand_help}")
        sys.exit(1)


def setup_frontend(
        package_main_name: str,
        one_line_description: str,
        version: str,
        help_info: str=None,
        subcommand_help: str = "Use 'lscmd' to list all valid subcommands."
):
    if os.environ.get('LOG_LEVEL') is None:
        logger_helper.set_level('INFO')
    lh.info(f'{one_line_description} ver. {version}')
    lh.info(f'Called by: {" ".join(sys.argv)}')
    parsed_args = _parse_args(sys.argv[1:])
    _act_on_args(
        parsed_args,
        package_main_name,
        version,
        help_info,
        subcommand_help
    )
