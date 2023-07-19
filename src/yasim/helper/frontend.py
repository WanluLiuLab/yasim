"""
frontend.py -- Utilities for other ``main`` functions.

.. versionadded:: 3.1.5
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
    help="Depth lower than this value would be this value.",
    nargs='?',
    type=float,
    action='store',
    default=0.01
))
FrontendOptSpecs.add(FrontendOptSpec(
    '--high_cutoff_ratio',
    required=False,
    help="Depth higher than `mu * high_cutoff_ratio` would be `mu * high_cutoff_ratio`",
    nargs='?',
    type=float,
    action='store',
    default=200
))
FrontendOptSpecs.add(FrontendOptSpec(
    '--preserve_intermediate_files',
    required=False,
    help="Do not remove intermediate files.",
    action='store_true',
))


def patch_frontend_argument_parser(parser: argparse.ArgumentParser, argname: str) -> argparse.ArgumentParser:
    """
    TODO docs
    
    .. versionadded:: 3.1.5
    """
    return FrontendOptSpecs.patch(parser, argname)
