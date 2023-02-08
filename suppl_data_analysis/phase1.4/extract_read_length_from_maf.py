import re
import sys
from typing import Tuple, Iterable

from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader

MafRecordType = Tuple[str, str, str, str]
"""Name1, Name2, Alignment1, Alignment2"""

maf_record_regex = re.compile(r"^s +(\S+) +([0-9]+) +([0-9]+) +([+-]) +([0-9]+) +(\S+)$")


def maf_parse(maf_path: str) -> Iterable[MafRecordType]:
    num_error = 0
    num_record = 0
    with get_tqdm_line_reader(maf_path) as reader:
        while True:
            line1 = reader.readline()
            if line1 == "":
                break
            if not line1.startswith("s"):
                continue
            else:
                line2 = reader.readline().rstrip("\n")
                line1 = line1.rstrip("\n")
            lm1 = maf_record_regex.match(line1)
            lm2 = maf_record_regex.match(line2)
            if lm1 is None or lm2 is None:
                num_error += 1
                continue
            g1 = lm1.groups()
            g2 = lm2.groups()
            yield g1[0], g2[0], g1[5], g2[5]
            num_record += 1
    print(f"Finished with {num_error} errors and {num_record} records", file=sys.stderr)


if __name__ == "__main__":
    all_qual = {"I":0, "D":0, "M":0, "S":0}
    print("\t".join(("ALIGNED_TRANSCRIPT_ID", "SIMULATED_TRANSCRIPT_ID", "READ_LENGTH")))
    for maf_record in maf_parse(sys.argv[1]):
        aligned_transcript_id = maf_record[0]
        simulated_transcript_id = maf_record[1].split(":")[0]
        read_length = len(maf_record[3].replace("-", ""))
        print("\t".join((
            aligned_transcript_id,
            simulated_transcript_id,
            str(read_length)
        )))
