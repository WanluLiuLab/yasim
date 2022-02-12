import argparse
import sys
from typing import List

from bioutils.datastructure import GeneView

from bioutils.datastructure.gene_typing import Transcript
from commonutils import ioctl

__version__ = 0.1

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
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))
    return parser.parse_args(args)

# FIXME
def get_with_defaults(transcript: Transcript, key: str) -> str:
    if key in transcript.data.attribute.keys():
        return repr(transcript.data.attribute[key])
    else:
        return "NA"


def main(args: List[str]):
    args = _parse_args(args)
    ans_gv = GeneView.from_file(args.gtf)
    with ioctl.get_writer(args.out) as writer:
        writer.write("\t".join(POSSIBLE_KEYS) + "\n")
        for transcript in ans_gv.transcripts.values():
            current_line = []
            for key in POSSIBLE_KEYS:
                current_line.append(get_with_defaults(transcript, key))
            writer.write("\t".join(current_line) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
