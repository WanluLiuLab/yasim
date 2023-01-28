import re
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
    print(f"Finished with {num_error} errors and {num_record} records")


if __name__ == "__main__":
    all_qual = {"I":0, "D":0, "M":0, "S":0}
    for maf_record in maf_parse("simulated.maf.gz"):
        # LAST Generated MAF         {'I': '6.05%', 'D': '4.13%', 'M': '70.81%', 'S': '19.01%'} 869316 lines
        # PBSIM Generated MAF        {'I': '7.71%', 'D': '5.46%', 'M': '85.98%', 'S': '0.84%'}  869316 lines
        # PBSIM Generated Log (Mean) {'I': '8.12%', 'D': '5.74%', 'M': '85.23%', 'S': '0.88%'}
        for b1, b2 in zip(maf_record[2], maf_record[3]):
            if b1 == "-":
                this_cigar = "I"
            elif b2 == "-":
                this_cigar= "D"
            elif b1 == b2:
                this_cigar= "M"
            else:
                this_cigar= "S"
            all_qual[this_cigar] += 1
    print("Events:", {k:str(round(v / sum(all_qual.values())*100, 2)) + "%" for k, v in all_qual.items()})

