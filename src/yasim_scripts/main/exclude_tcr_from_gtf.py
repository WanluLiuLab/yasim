from typing import List

from labw_utils.bioutils.parser.feature import GtfIterator, GtfWriter


def main(args: List[str]):
    for arg in args:
        with GtfWriter(arg + ".tcr_filtered.gtf") as writer:
            for record in GtfIterator(arg, show_tqdm=True):
                if not record.attribute.get("gene_name", "TRBV").startswith("TR"):
                    writer.write_feature(record)
