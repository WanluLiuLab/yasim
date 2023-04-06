"""
frontend.py -- Utilities for other ``main`` functions.
"""

__all__ = (
    "patch_frontend_argument_parser",
)

import argparse
from typing import Dict, Callable


def _add_f(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        '-f',
        '--fasta',
        required=True,
        help="Path to input reference genome sequence, in FASTA format",
        nargs='?',
        type=str,
        action='store'
    )
    return parser


def _add_g(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        '-g', '--gtf',
        required=True,
        help="Path to input genomic annotation in GTF format. Can be compressed.",
        nargs='?',
        type=str,
        action='store'
    )
    return parser


def _add_d(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        '-d', '--depth',
        required=True,
        help="Path to input Isoform-Level Depth TSV. Can be compressed.",
        nargs='?',
        type=str,
        action='store'
    )
    return parser


def _add_b(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        '-b', '--barcodes',
        required=True,
        help="Path to input barcode TXT. Can be compressed.",
        nargs='?',
        type=str,
        action='store'
    )
    return parser


def _add_lc(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser.add_argument(
        '--low_cutoff',
        required=False,
        help="Depth lower than this value would be 0.",
        nargs='?',
        type=float,
        action='store',
        default=0.01
    )
    return parser


_PARSER_DICT: Dict[
    str,
    Callable[[argparse.ArgumentParser], argparse.ArgumentParser]
] = {
    "-f": _add_f,
    "-g": _add_g,
    "-d": _add_d,
    "--low_cutoff": _add_lc,
    "-b": _add_b
}


def patch_frontend_argument_parser(parser: argparse.ArgumentParser, argname: str) -> argparse.ArgumentParser:
    return _PARSER_DICT[argname](parser)
