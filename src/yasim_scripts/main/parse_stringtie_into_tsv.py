"""Parse stringtie GTF output into a TSV file suitable for downstream analysis using R"""

import argparse
from typing import List

from bioutils.io.feature import GtfIterator
from bioutils.typing.feature import GtfRecord
from commonutils.io.safe_io import get_writer

POSSIBLE_KEYS = (
    "gene_id",
    "transcript_id",
    "reference_id",
    "ref_gene_id",
    "ref_gene_name",
    "cov",
    "FPKM",
    "TPM"
)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', required=True, help="Stringtie Output GTF", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output TSV", nargs='?',
                        type=str, action='store')
    return parser.parse_args(args)


def get_with_defaults(transcript_gtf_record: GtfRecord, key: str) -> str:
    if key in transcript_gtf_record.attribute.keys():
        return repr(transcript_gtf_record.attribute[key])
    else:
        return "NA"


def main(args: List[str]):
    args = _parse_args(args)
    ans_gv = GtfIterator(args.gtf)
    with get_writer(args.out) as writer:
        writer.write("\t".join(POSSIBLE_KEYS) + "\n")
        for gtf_record in ans_gv:
            if gtf_record.feature == "transcript":
                current_line = []
                for key in POSSIBLE_KEYS:
                    current_line.append(get_with_defaults(gtf_record, key))
                writer.write("\t".join(current_line) + "\n")
