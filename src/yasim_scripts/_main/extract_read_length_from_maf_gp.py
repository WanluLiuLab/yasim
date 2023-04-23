"""
extract_read_length_from_maf_gp.py -- Extraction of Read Length from MAF, General-Purposed

This script can be used to extract read length of all transcript ID from transcriptomically-aligned TGS RNA-Seq MAF
for assessing read completeness.

Synopsis: python -m yasim_scrips extract_quality_from_maf [MAF1] [[MAF2]...]

Arguments:
    [MAF1] [[MAF2]...] path to MAF files produced by PBSIM3 or LAST aligner.
"""
__all__ = (
    "main",
)

from labw_utils.typing_importer import List

from yasim_scripts.helper.maf_parser import maf_parse


def main(args: List[str]):
    for arg in args:
        print("\t".join(("ALIGNED_TRANSCRIPT_ID", "READ_LENGTH")))
        for maf_record in maf_parse(arg):
            aligned_transcript_id = maf_record[0]
            read_length = len(maf_record[3].replace("-", ""))
            print("\t".join((
                aligned_transcript_id,
                str(read_length)
            )))
