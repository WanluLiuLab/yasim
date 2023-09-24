"""

"""
import argparse
import json
from collections import defaultdict

import pysam
import pandas as pd

from labw_utils.commonutils.appender import TableAppenderConfig, load_table_appender_class, BaseTableAppender
from labw_utils.commonutils.lwio import get_reader
from labw_utils.commonutils.stdlib_helper.argparse_helper import ArgumentParserWithEnhancedFormatHelp
from labw_utils.typing_importer import Dict, Set, List
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from yasim.helper.translation_instruction import TranslationInstruction, SimpleTE
from yasim.helper.rmsk_parser import RMSKGffIterator

ReadIDTENameMap = Dict[str, Set[str]]


def convert_hmmer_to_tes(src_hmmer_path: str) -> ReadIDTENameMap:
    retd = defaultdict(lambda: set())
    with get_reader(src_hmmer_path) as r:
        for line in r:
            if line.startswith("#"):
                continue
            try:
                ls = list(filter(lambda x: bool(x), line.strip().split(" ")))
                if float(ls[12]) < 1E-5:
                    retd[ls[0]].add(ls[1].split("#")[0])
            except Exception:
                print(f"ERR: {line}")
    return retd

def convert_aln_bam_to_tes(src_aln_bam_path: str, is_loci: bool) -> ReadIDTENameMap:
    retd = defaultdict(lambda: set())
    with pysam.AlignmentFile(src_aln_bam_path) as af:
        for record in af.fetch():
            te_name = record.reference_name
            if is_loci:
                te_name = "_".join(te_name.split("_")[:-1])
            retd[record.query_name].add(te_name)
    return retd


def convert_aln_blast6_to_tes(src_aln_blast6_path: str, is_loci: bool) -> ReadIDTENameMap:
    df = (
        pd.read_csv(
            src_aln_blast6_path,
            sep="\t",
            names=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "mismatch",
                "gapopen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "evalue",
                "bitscore",
            ],
            engine="c",
            comment="#",
            usecols=["qseqid", "sseqid", "evalue", "length"],
        )
        .query("evalue < 1E-5")
        .query("length > 20")
    )
    if is_loci:
        df["sseqid"] = df["sseqid"].apply(lambda s: "_".join(s.split("_")[:-1]))
    retd = defaultdict(lambda: set())
    for tup in df.itertuples(index=False):
        retd[tup.qseqid].add(tup.sseqid)
    return retd


def convert_rmsk_gff_to_tes(src_aln_rmsk_gff_path: str) -> ReadIDTENameMap:
    retd = defaultdict(lambda: set())
    with RMSKGffIterator(src_aln_rmsk_gff_path) as gffi:
        for record in gffi:
            repeat_name = record.attribute_get("repeat_name")
            if repeat_name.endswith(")n") or repeat_name.endswith("-rich"):
                continue
            if record.end0b - record.start0b + 1 < 20:
                continue
            retd[record.seqname].add(record.attribute_get("repeat_name"))
    return retd


def convert_fc_assignment_to_tes(src_fc_assignment_path: str) -> ReadIDTENameMap:
    df = pd.read_csv(
        src_fc_assignment_path,
        sep="\t",
        names=[
            "qseqid",
            "status",
            "unknown",
            "sseqid",
        ],
        dtype={
            "qseqid": str,
            "status": str,
            "unknown": str,
            "sseqid": str,
        },
        engine="c",
        comment="#",
    ).query('status == "Assigned"')
    df["sseqid"] = df["sseqid"].apply(lambda s: "_".join(s.split("_")[:-1]))
    retd = defaultdict(lambda: set())
    for tup in df.itertuples(index=False):
        retd[tup.qseqid].add(tup.sseqid)
    return retd


def convert_translation_instruction_to_tes(src_translation_instruction_path: str) -> ReadIDTENameMap:
    ti = TranslationInstruction.from_json(src_translation_instruction_path)
    return {
        k: set(map(lambda _x: _x.src_te_name, filter(lambda x: isinstance(x, SimpleTE), v.l)))
        for k, v in ti.transcripts.items()
    }


def comp_tes_tes(src_gt_tes_file_path: str, src_query_tes_file_path: str, dst_path: str):
    with get_reader(src_gt_tes_file_path, is_binary=False) as r:
        gt_tes = json.load(r)
    with get_reader(src_query_tes_file_path, is_binary=False) as r:
        query_tes = json.load(r)

    appender: BaseTableAppender = load_table_appender_class("TSVTableAppender")(
        dst_path,
        ["seq_id", "te_name", "status"],
        tac=TableAppenderConfig(),
    )
    for gt_id in tqdm(gt_tes.keys()):
        query_te_names = set(query_tes.get(gt_id, set()))
        gt_te_names = set(gt_tes[gt_id])
        all_te_names = query_te_names.union(gt_te_names)
        for te_name in all_te_names:
            if te_name in gt_te_names and te_name in query_te_names:
                appender.append([gt_id, te_name, "TP"])
            elif te_name in gt_te_names:
                appender.append([gt_id, te_name, "FN"])
            elif te_name in query_te_names:
                appender.append([gt_id, te_name, "FP"])
    appender.flush()
    appender.close()


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim insert_transposon",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument("--src_gt_tes_file_path", type=str, required=True)
    parser.add_argument("--src_query_tes_file_path", type=str, required=True)
    parser.add_argument("--dst_path", type=str, required=True)
    return parser


def main(args: List[str]):
    argv = create_parser().parse_args(args)
    comp_tes_tes(
        src_gt_tes_file_path=argv.src_gt_tes_file_path,
        src_query_tes_file_path=argv.src_query_tes_file_path,
        dst_path=argv.dst_path,
    )
