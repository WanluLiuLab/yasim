"""
rearrange_tcr.py -- Generate ground-truth TCR Contigs
"""
__all__ = (
    "main",
    "create_parser"
)

import argparse
import json
import os
from typing import List, Dict

from labw_utils.commonutils.appender import load_table_appender_class, TableAppenderConfig
from labw_utils.commonutils.io.safe_io import get_reader, get_writer
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from yasim.helper.frontend import patch_frontend_argument_parser
from yasim.helper.tcr import TCRTranslationTableType, Cdr3InsertionTable, TCell, GenerationFailure
from yasim_sctcr._main import get_sample_data_path

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_sctcr rearrange_tcr",
        description=__doc__.splitlines()[1]
    )
    parser.add_argument(
        '--tcr_genelist_path',
        required=True,
        help=f"TCR Gene List JSON. The IMGT version for human is {get_sample_data_path('tcr_gene_list.min.json')}",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '--tcr_cache_path',
        required=True,
        help="TCR Cache JSON generated by `generate_tcr_cache`",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '--cdr3_deletion_table_path',
        required=True,
        help=f"TCR CDR3 Deletion Table JSON. See examples in {get_sample_data_path('cdr3_deletion_table.min.json')}",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '--cdr3_insertion_table_path',
        required=True,
        help=f"TCR CDR3 Insertion Table JSON. See examples in {get_sample_data_path('cdr3_insertion_table.min.json')}",
        nargs='?',
        type=str,
        action='store'
    )
    parser = patch_frontend_argument_parser(parser, "-b")
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Output basename",
        nargs='?',
        type=str,
        action='store'
    )
    return parser


def main(args: List[str]):
    args = create_parser().parse_args(args)
    rearrange_tcr(
        output_base_path=args.out,
        tcr_genelist_path=args.tcr_genelist_path,
        tcr_cache_path=args.tcr_cache_path,
        cdr3_deletion_table_path=args.cdr3_deletion_table_path,
        cdr3_insertion_table_path=args.cdr3_insertion_table_path,
        barcode_path=args.barcodes
    )


def rearrange_tcr(
        barcode_path: str,
        tcr_cache_path: str,
        tcr_genelist_path: str,
        output_base_path: str,
        cdr3_deletion_table_path: str,
        cdr3_insertion_table_path: str
):
    n_failure = 0
    with get_reader(tcr_genelist_path) as reader:
        tcr_genelist = json.load(reader)
    with get_reader(cdr3_deletion_table_path) as reader:
        cdr3_deletion_table = json.load(reader)
    with get_reader(tcr_cache_path) as reader:
        tcr_cache: Dict[str, TCRTranslationTableType] = json.load(reader)
    cdr3_insertion_table = Cdr3InsertionTable(cdr3_insertion_table_path)
    cdr3_deletion_table = {
        k: {
            int(_k): _v
            for _k, _v in v.items()
        }
        for k, v in cdr3_deletion_table.items()
    }
    with get_writer(output_base_path + ".nt.fa") as nt_fasta_writer, \
            get_writer(output_base_path + ".aa.fa") as aa_fasta_writer:
        with load_table_appender_class("TSVTableAppender")(
                filename=output_base_path + ".stats",
                header=[
                    "UUID",
                    "TRAV",
                    "TRAJ",
                    "TRBV",
                    "TRBJ",
                    "ACDR3_AA",
                    "BCDR3_AA",
                    "ACDR3_NT",
                    "BCDR3_NT",
                    "ALPHA_AA",
                    "BETA_AA",
                    "ALPHA_NT",
                    "BETA_NT"
                ],
                tac=TableAppenderConfig(buffer_size=1024)
        ) as appender:
            for barcode in get_tqdm_line_reader(barcode_path):
                while True:
                    try:
                        cell = TCell.from_gene_names(
                            tcr_genelist=tcr_genelist,
                            cdr3_deletion_table=cdr3_deletion_table,
                            cdr3_insertion_table=cdr3_insertion_table,
                            tcr_cache=tcr_cache,
                            barcode=barcode
                        )
                    except GenerationFailure:
                        n_failure += 1
                        continue
                    else:
                        break
                nt_fasta_writer.write(
                    cell.to_nt_fasta_record() + "\n"
                )
                aa_fasta_writer.write(
                    cell.to_aa_fasta_record() + "\n"
                )
                appender.append([
                    cell.cell_uuid,
                    *cell.alpha_names,
                    *cell.beta_names,
                    *cell.cdr3_aa,
                    *cell.cdr3_nt,
                    cell.alpha_aa,
                    cell.beta_aa,
                    cell.alpha_nt,
                    cell.beta_nt
                ])
                with get_writer(os.path.join(output_base_path + ".json.d", cell.cell_uuid + ".json")) as writer:
                    json.dump(cell.to_dict(), writer)
    _lh.info("Finished with %d failures", n_failure)
