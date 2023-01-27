import re
from typing import Tuple, Iterable

from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader

MafRecordType = Tuple[str, str, str, str]
"""Name1, Name2, Alignment1, Alignment2"""

maf_record_regex = re.compile(r"^s +(\S+) +([0-9]+) +([0-9]+) +([+-]) +([0-9]+) +(\S+)$")


def maf_parse(maf_path: str) -> Iterable[MafRecordType]:
    tmpl = []
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
    all_fake_cigar = ""
    for maf_record in maf_parse("/home/yuzj/Documents/yasim/suppl_data_analysis/phase1.3/simulated.maf"):
        # {'I': 9.99, 'D': 3.79, 'M': 80.2, 'S': 6.01}
        # {'I': 9.54, 'D': 3.68, 'M': 80.88, 'S': 5.9}
        fake_cigar = ""
        for b1, b2 in zip(maf_record[2], maf_record[3]):
            if b1 == "-":
                fake_cigar += "I"
            elif b2 == "-":
                fake_cigar += "D"
            elif b1 == b2:
                fake_cigar += "M"
            else:
                fake_cigar += "S"
        # print(maf_record[2] + "\n" + maf_record[3] + "\n" + fake_cigar + "\n")
        all_fake_cigar += fake_cigar
    all_qual = {
        k: all_fake_cigar.count(k)
        for k in "IDMS"
    }
    print("Events:", {k:str(round(v / len(all_fake_cigar)*100, 2)) + "%" for k, v in all_qual.items()})

