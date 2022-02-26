"""

To parse some file looks like:

SIRV1	11404	11606	+	gene_id SIRV1; transcript_id SIRV108; exon_assignment SIRV108_2; count=14

into:

SIRV1	NanoAsPipe	exon	11404	11606	.	+	.	gene_id SIRV1; transcript_id SIRV108; exon_assignment SIRV108_2; count=14
"""

import argparse
from typing import List

__version__ = 0.1


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, help="NanoAsPipe input file", nargs='?',
                        type=str, action='store')
    parser.add_argument('-o', '--out', required=True, help="Output TSV", nargs='?',
                        type=str, action='store')
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)

