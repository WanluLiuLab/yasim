import glob
import os

from bioutils.io.feature import Gff3Iterator


if __name__ == "__main__":
    for data_filename in glob.glob(os.path.join(os.path.dirname(__file__),"data","*.gff3")):
        gffi = Gff3Iterator(data_filename)
        for record in gffi:
            if record.feature == "gene":
                print("\t".join((
                    record.attribute['GeneName'],
                    record.source
                    )))
