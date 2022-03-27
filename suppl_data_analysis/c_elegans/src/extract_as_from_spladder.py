import glob
import os

from bioutils.io.feature import Gff3Iterator
from commonutils.io.safe_io import get_writer


if __name__ == "__main__":
    with get_writer("all_as.tsv") as writer:
        for data_filename in glob.glob(os.path.join(os.path.dirname(__file__), "..","as_data", "*.gff3")):
            gffi = Gff3Iterator(data_filename)
            for record in gffi:
                if record.feature == "gene":
                    writer.write("\t".join((
                        record.attribute['GeneName'],
                        record.source
                        ))+"\n")

