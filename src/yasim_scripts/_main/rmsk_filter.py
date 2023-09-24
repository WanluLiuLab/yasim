"""
Filter and convert RepeatMasker GFF to GTF
"""
import argparse

from labw_utils.bioutils.parser.gtf import GtfIteratorWriter
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.typing_importer import List
from yasim.helper import rmsk_parser


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_scripts rmsk_filter",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument(
        "-i",
        "--src_rmsk_gff_path",
        type=str,
        required=True,
        help="File path of source RepeatMasker GFF",
    )
    parser.add_argument(
        "-o",
        "--dst_rmsk_gtf_path",
        type=str,
        required=True,
        help="File path of destination TE GFF",
    )
    parser.add_argument("--no_filter_simple_repeats", action="store_false")
    parser.add_argument(
        "--filter_short_repeats_threshold",
        type=int,
        default=20,
        required=False,
    )
    return parser


def filter_rmsk(
    src_rmsk_gff_path: str,
    dst_rmsk_gtf_path: str,
    filter_simple_repeats: bool,
    filter_short_repeats_threshold: str,
):
    with GtfIteratorWriter(dst_rmsk_gtf_path) as gtfw:
        with rmsk_parser.RMSKGffIterator(src_rmsk_gff_path, True) as rmski:
            for i, feature in enumerate(rmski):
                repeat_name = feature.attribute_get("repeat_name")
                if filter_simple_repeats and (repeat_name.endswith(")n") or repeat_name.endswith("-rich")):
                    continue
                if feature.end0b - feature.start0b + 1 < filter_short_repeats_threshold:
                    continue
                gtfw.write(feature.update_attribute(gene_id=repeat_name, transcript_id=repeat_name + f"_{i}"))


def main(args: List[str]):
    argv = create_parser().parse_args(args)
    filter_rmsk(
        src_rmsk_gff_path=argv.src_rmsk_gff_path,
        dst_rmsk_gtf_path=argv.dst_rmsk_gtf_path,
        filter_simple_repeats=not argv.filter_simple_repeats,
        filter_short_repeats_threshold=argv.filter_short_repeats_threshold,
    )
