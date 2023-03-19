import argparse
import random
from typing import List, Set

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_writer
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

logger = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n',
        '--num_cells',
        required=True,
        help="Number of Cells",
        nargs='?',
        type=int,
        action='store'
    )
    parser.add_argument(
        '-l', '--length',
        required=False,
        default=17,
        help="Length of barcode",
        nargs='?',
        type=int,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Output barcode TXT",
        nargs='?',
        type=str,
        action='store'
    )
    return parser.parse_args(args)


class BarcodeGenerator:
    _generated_barcodes: Set[str] = set()

    @staticmethod
    def generate(barcode_length: int) -> str:
        while True:
            new_barcode = "".join(map(lambda _: random.choice("ATCG"), range(barcode_length)))
            if new_barcode not in BarcodeGenerator._generated_barcodes:
                BarcodeGenerator._generated_barcodes.add(new_barcode)
                return new_barcode


def main(args: List[str]):
    args = _parse_args(args)
    with get_writer(args.out) as writer:
        for _ in tqdm(range(args.num_cells), desc="Generating"):
            writer.write(BarcodeGenerator.generate(args.length) + "\n")
