import itertools
import json
import os
from typing import List, Tuple

import numpy as np

from labw_utils.bioutils.algorithm.alignment import SmithWatermanAligner
from labw_utils.bioutils.algorithm.sequence import is_valid_chrname, translate_cdna
from labw_utils.bioutils.datastructure.fasta_view import FastaViewFactory
from labw_utils.bioutils.datastructure.gene_view import GeneViewFactory
from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_reader, get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

lh = get_logger(__name__)

TCRTranslationTableType = List[Tuple[str, str, str, str]]
"""
[(REAL_NT, TRANS_AA, REAL_AA, STATUS)]
"""


def align(real_nt_seq: str, real_aa_seq: str) -> TCRTranslationTableType:
    swas: List[SmithWatermanAligner] = []
    translation_table = []
    for i in (0, 1, 2):
        seq = real_nt_seq[i:]
        seq = seq[: len(seq) - len(seq) % 3]
        swa = SmithWatermanAligner(
            translate_cdna(seq),
            real_aa_seq,
            is_global=False,
            mismatch_score=-6,
            indel_score=-6
        )
        swas.append(swa)
    mn = np.argmax(list(swa.score for swa in swas))
    backtrack = swas[mn].get_backtrack()[0].splitlines()[1:]
    """
    Should be like:
    
    MRLVARVTVFLTFGTIIDAKTTQPTSMDCAEGRAANLPCNHSTISGNEYVYWYRQIHSQGPQYIIHGLKNNETNEMASLIITEDRKSSTLILPHATLRDTAVYYCIVRV
    DDDDDDDDDDDDDDDDD=======M====================================================================================
    -----------------DAKTTQPPSMDCAEGRAANLPCNHSTISGNEYVYWYRQIHSQGPQYIIHGLKNNETNEMASLIITEDRKSSTLILPHATLRDTAVYYCIVRV
    """
    nt_p = 0
    if mn != 0:
        nt_p = mn
        translation_table.append((real_nt_seq[0:nt_p], "-", "-", "X"))
    for trans_aa, status, real_aa in zip(*backtrack):
        if status == "=" or status == "M":
            translation_table.append((real_nt_seq[nt_p:nt_p + 3], trans_aa, real_aa, status))
            nt_p += 3
        elif status == "I":
            translation_table.append(("-", "-", real_aa, status))
        elif status == "D":
            translation_table.append((real_nt_seq[nt_p:nt_p + 3], trans_aa, "-", status))
            nt_p += 3
    if nt_p < len(real_nt_seq):
        translation_table.append((real_nt_seq[nt_p:], "-", "-", "X"))
    return translation_table


def create_tcr_cache(
        ref_fa_path: str,
        ref_gtf_path: str,
        tcr_genelist_path: str,
        tcr_aa_table_path: str,
        tcr_cache_path: str
):
    with get_reader(tcr_genelist_path) as reader:
        tcr_genelist = json.load(reader)
    with get_reader(tcr_aa_table_path) as reader:
        tcr_aa_table = json.load(reader)
    ref_fasta_view = FastaViewFactory(
        ref_fa_path,
        show_tqdm=True,
        read_into_memory=False
    )
    ref_gene_view = GeneViewFactory.from_file(ref_gtf_path)
    tcrs = {}
    for gene_name in tqdm(
            list(itertools.chain(*tcr_genelist.values())),
            desc="Generating TCR Cache"
    ):
        try:
            gene = ref_gene_view.get_gene(gene_name)
        except KeyError:
            lh.warning("Gene %s not found!", gene_name)
            continue
        valid_transcript = []
        for transcript in gene.iter_transcripts():
            if not is_valid_chrname(transcript.seqname):
                continue
            valid_transcript.append(transcript)
        if len(valid_transcript) != 1:
            lh.warning(
                "Gene %s have %d != 1 valid transcript (%s)!",
                gene_name,
                len(valid_transcript),
                str(list(map(lambda x: x.transcript_id, valid_transcript)))
            )
            continue
        real_nt_seq = valid_transcript[0].cdna_sequence(ref_fasta_view.sequence)
        if real_nt_seq == "":
            lh.warning(
                "Gene %s have no NT sequence!",
                gene_name
            )
            continue
        try:
            real_aa_seq = tcr_aa_table["human"][gene_name]
        except KeyError:
            lh.warning(
                "Gene %s have no AA sequence!",
                gene_name
            )
            continue
        tcrs[gene_name] = align(real_nt_seq=real_nt_seq, real_aa_seq=real_aa_seq)
    with get_writer(tcr_cache_path) as writer:
        json.dump(tcrs, writer)

    with get_writer(tcr_cache_path + ".fa") as writer:
        for tcr_name, tcr_tt in tcrs.items():
            writer.write(f">{tcr_name}\n{''.join(list(zip(*tcr_tt))[0])}\n")


if __name__ == "__main__":
    basename = "."
    create_tcr_cache(
        ref_fa_path=os.path.join(basename, "hg38.fa"),
        ref_gtf_path=os.path.join(basename, "hg38.ncbiRefSeq_chr7_14.gtf"),
        tcr_genelist_path=os.path.join(basename, "tcr_genelist.json"),
        tcr_aa_table_path=os.path.join(basename, "IMGT_Protein_Display.json"),
        tcr_cache_path=os.path.join(basename, "tcr_cache.json")
    )
