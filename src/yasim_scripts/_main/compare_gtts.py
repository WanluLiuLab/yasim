"""

"""
from collections import defaultdict

import pysam
import pandas as pd

from labw_utils.typing_importer import Dict, Set, List
from yasim.helper.translation_instruction import TranslationInstruction, SimpleTE
from yasim.helper.rmsk_parser import RMSKGffIterator


ReadIDTENameMap: Dict[str, Set[str]]


def convert_aln_bam_to_tes(src_aln_bam_path: str, is_loci: bool) -> ReadIDTENameMap:
    retd = defaultdict(lambda: set())
    with pysam.AlignmentFile(src_aln_bam_path) as af:
        for record in af.fetch():
            ...
    return retd


def convert_aln_blast6_to_tes(src_aln_blast6_path: str, is_loci: bool) -> ReadIDTENameMap:
    df = pd.read_csv(
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
        usecols=["qseqid", "sseqid", "evalue"],
    ).query("evalue < 1E-5")
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
        engine="c",
        comment="#",
    ).query('sseqid != "NA"')
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


def main(args: List[str]):
    pass
