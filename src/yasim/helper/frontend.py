"""
frontend.py -- Utilities for other ``main`` functions.
"""

__all__ = (
    "patch_frontend_argument_parser",
)

import argparse

from labw_utils.bioutils.comm_frontend_opts import FrontendOptSpecs, FrontendOptSpec

FrontendOptSpecs.add(FrontendOptSpec(
    '-d', '--depth',
    required=True,
    help="Path to input Isoform-Level Depth TSV. Can be compressed.",
    nargs='?',
    type=str,
    action='store'
))

FrontendOptSpecs.add(FrontendOptSpec(
    '-b', '--barcodes',
    required=True,
    help="Path to input barcode TXT. Can be compressed.",
    nargs='?',
    type=str,
    action='store'
))
FrontendOptSpecs.add(FrontendOptSpec(
    '--low_cutoff',
    required=False,
    help="Depth lower than this value would be 0.",
    nargs='?',
    type=float,
    action='store',
    default=0.01
))


def patch_frontend_argument_parser(parser: argparse.ArgumentParser, argname: str) -> argparse.ArgumentParser:
    return FrontendOptSpecs.patch(argname, parser)
