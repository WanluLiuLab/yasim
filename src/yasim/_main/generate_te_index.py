"""
generate_te_index.py -- Generate Transposon Indices from DFam HDF5 file

.. versionadded:: 3.2.0
"""

__all__ = ("create_parser", "main")

import argparse

from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.typing_importer import List
from yasim.helper.transposon import TransposonDatabase


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim generate_te_index",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument(
        "-i",
        "--src_dfam_h5_file_path",
        required=True,
        help="Input Dfam HDF5 file path. Get one from <https://dfam.org/releases/current/families/>",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "-o",
        "--dst_te_idx_file_path",
        required=True,
        help="Output TE index path",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "--dst_consensus_fa_path",
        required=False,
        default=None,
        help="Output TE index path",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "--dst_hmm_path",
        required=False,
        default=None,
        help="Output TE index path",
        nargs="?",
        type=str,
        action="store",
    )
    parser.add_argument(
        "-t",
        "--txid",
        required=False,
        help="NCBI taxon ID to subset.",
        nargs="?",
        type=int,
        action="store",
        default=9606,
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    TransposonDatabase.convert_dfam_hdf5(
        src_dfam_hdf5_file_path=args.src_dfam_h5_file_path,
        dst_index_file_path=args.dst_te_idx_file_path,
        dst_consensus_fa_path=args.dst_consensus_fa_path,
        dst_hmm_path=args.dst_hmm_path,
        with_tqdm=True,
        required_txid=args.txid,
    )
