"""
insert_transposon.py -- Generate 

.. versionadded:: 3.2.0
"""

__all__ = (
    "main",
    "create_parser"
)

import argparse
import json
import random
from labw_utils.commonutils.lwio.safe_io import get_writer

from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.bioutils.parser.fasta import FastaIterator, FastaWriter
from labw_utils.bioutils.record.fasta import FastaRecord
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List

from yasim.helper.transposon import TransposonDatabase

rdg = random.SystemRandom()

def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim insert_transposon",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument(
        "--te_dbi",
        type=str,
        help="Path to TE database index.",
        required=True
    )
    parser.add_argument(
        "--src_fasta_path",
        type=str,
        help="Path to FASTA where TEs would be inserted into.",
        required=True   
    )
    parser.add_argument(
        "--dst_fasta_path",
        type=str,
        help="Path to FASTA where result file should be written to.",
        required=True   
    )
    parser.add_argument(
        "--dst_frag_spec_path",
        type=str,
        help="Path to generated ground-truth fragment specifications",
        required=True   
    )
    parser.add_argument(
        "--percent_reads_add_te",
        type=float,
        help="Percentile of reads with TEs",
        required=False,
        default= 40.0
    )
    return parser




def main(args: List[str]):
    def whether_to_generate_te() -> bool:
        return random.choices(
            [True, False],
            weights=(percent_reads_add_te, 100 - percent_reads_add_te),
            k = 1
        )[0]

    _lh = get_logger(__name__)
    argv = create_parser().parse_args(args)
    if not argv.percent_reads_add_te < 100:
        _lh.error(f"percent_reads_add_te should be between 0 and 100, currently {argv.percent_reads_add_te}")
        return 1
    percent_reads_add_te = argv.percent_reads_add_te
    tedb = TransposonDatabase.load(argv.te_dbi, with_tqdm=True)
    fragspec = {}
    with FastaIterator(argv.src_fasta_path, with_tqdm=True) as fai:
        with FastaWriter(argv.dst_fasta_path) as faw:
            for record in fai:
                fragspec[record.seq_id] = {}
                if whether_to_generate_te():
                    seq = record.sequence
                    ...
                    faw.write(FastaRecord(seq_id=record.seq_id, sequence=seq))
                else:
                    faw.write(record=record)
    with get_writer(argv.dst_frag_spec_path, is_binary = False) as w:
        json.dump(fragspec, w, indent=4)



