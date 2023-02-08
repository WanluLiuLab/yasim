import re
import sys
from typing import Tuple, Iterable, List

from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader

from yasim_scripts.main.extract_quality_from_maf import maf_parse

MafRecordType = Tuple[str, str, str, str]
"""Name1, Name2, Alignment1, Alignment2"""

maf_record_regex = re.compile(r"^s +(\S+) +([0-9]+) +([0-9]+) +([+-]) +([0-9]+) +(\S+)$")


def main(args:List[str]):
    print("\t".join(("ALIGNED_TRANSCRIPT_ID", "SIMULATED_TRANSCRIPT_ID", "READ_LENGTH")))
    for maf_record in maf_parse(args[0]):
        aligned_transcript_id = maf_record[0]
        simulated_transcript_id = maf_record[1].split(":")[0]
        read_length = len(maf_record[3].replace("-", ""))
        print("\t".join((
            aligned_transcript_id,
            simulated_transcript_id,
            str(read_length)
        )))
