import argparse
import sys
from typing import List

from bioutils.datastructure.gene import GeneView

from commonutils.tqdm_importer import tqdm

__version__ = 0.1


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gtf', required=True, help="Stringtie Output GTF", nargs='?',
                        type=str, action='store')
    parser.add_argument('-r', '--ref_gtf', required=True, help="Reference GTF", nargs='?',
                        type=str, action='store')
    parser.add_argument('-t', '--ground_truth', required=False, help="Ground Truth GTF", nargs='?',
                        type=str, action='store')
    parser.add_argument('-v', '--version', help="Print version information", action='version',
                        version='%(prog)s ' + str(__version__))

    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    ans_gv = GeneView(args.gtf)
    ref_gv = GeneView(args.ref_gtf)
    if args.ground_truth:
        truth_gv = GeneView(args.ground_truth)
        number_of_transcripts_in_ground_truth = len(truth_gv.transcripts)
    else:
        truth_gv = None
        number_of_transcripts_in_ground_truth = 0
    number_of_transcripts_in_answer = len(ans_gv.transcripts)
    number_of_transcripts_in_reference = len(ref_gv.transcripts)
    number_of_transcripts_match_reference = 0
    number_of_transcripts_match_ground_truth = 0

    for transcript in tqdm(iterable=ans_gv.transcripts.values(), desc="Iterating"):
        if "reference_id" in transcript.data.attribute.keys():
            number_of_transcripts_match_reference += 1
        if truth_gv is not None:
            for transcript_truth in truth_gv.transcripts.values():
                if transcript_truth == transcript:
                    number_of_transcripts_match_ground_truth += 1
                    break
    print(f"number_of_transcripts_in_answer: {number_of_transcripts_in_answer}")
    print(f"number_of_transcripts_in_reference: {number_of_transcripts_in_reference}")
    print(f"number_of_transcripts_in_ground_truth: {number_of_transcripts_in_ground_truth}")

    ref_match_rate = round(number_of_transcripts_match_reference / number_of_transcripts_in_answer * 100, 2)
    print(
        f"number_of_transcripts_match_reference: {number_of_transcripts_match_reference}={ref_match_rate}%")

    truth_match_rate = round(number_of_transcripts_match_ground_truth / number_of_transcripts_in_answer * 100, 2)
    if number_of_transcripts_in_ground_truth != 0:
        print(
            f"number_of_transcripts_match_ground_truth: {number_of_transcripts_match_ground_truth}={truth_match_rate}%")


if __name__ == "__main__":
    main(sys.argv[1:])
