import re
import sys
from typing import Tuple, Iterable, List

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
    print(
        f"Finished with {num_error} errors and {num_record} records",
        file=sys.stderr
    )


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
        retl = [sys.argv[1]]
        retl.extend(map(str, all_qual.values()))
        print("\t".join(retl))
