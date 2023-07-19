"""
extract_quality_from_maf.py -- Extract per-base alignment status from MAF files.

This script extracts per-base alignment status from MAF files into a TSV file of one line,
which have 4 fields for Insertion, Deletion, Match and Substitution.

Synopsis: python -m yasim_scrips extract_quality_from_maf [MAF1] [[MAF2]...]

Arguments:
    [MAF1] [[MAF2]...] path to MAF files produced by PBSIM3 or LAST aligner.

.. versionadded:: 3.1.5
"""

__all__ = (
    "main",
)

from labw_utils.typing_importer import List

from yasim_scripts.helper.maf_parser import maf_parse


def main(args: List[str]):
    for arg in args:
        all_qual = {"I": 0, "D": 0, "M": 0, "S": 0}
        for maf_record in maf_parse(arg):
            for b1, b2 in zip(maf_record[2].upper(), maf_record[3].upper()):
                if b1 == "-":
                    this_cigar = "I"
                elif b2 == "-":
                    this_cigar = "D"
                elif b1 == "N" or b2 == "N" or b1 == b2:
                    this_cigar = "M"
                else:
                    this_cigar = "S"
                all_qual[this_cigar] += 1
        retl = [arg]
        retl.extend(map(str, all_qual.values()))
        print("\t".join(retl))
