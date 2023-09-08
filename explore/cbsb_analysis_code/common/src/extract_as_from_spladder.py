import argparse
import glob
import os

from labw_utils.bioutils.parser.feature import Gff3Iterator

from labw_utils.commonutils.lwio.safe_io import get_writer

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--spladder_dir", required=True)
    args = parser.parse_args()
    spladder_dir = args.spladder_dir
    with get_writer(os.path.join(spladder_dir, "all_as.tsv")) as writer:
        for data_filename in glob.glob(os.path.join(spladder_dir, "*.gff3")):
            gffi = Gff3Iterator(data_filename)
            for record in gffi:
                if record.feature == "gene":
                    writer.write("\t".join((record.attribute["GeneName"], record.source)) + "\n")
