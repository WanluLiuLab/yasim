import random

import pandas as pd

from labw_utils.bioutils.parser.feature import GtfIterator, GtfWriter
from labw_utils.bioutils.algorithm.sequence import is_valid_chrname

if __name__ == '__main__':
    df = pd.read_table("ncbi_dataset.tsv")
    pcg_df = (
        df.
        query("`Gene Type` == 'PROTEIN_CODING'").
        query("`Gene Group Method` == 'NCBI Ortholog'")
    )
    print(
        f"Filtered {len(pcg_df)} genes out of {len(df)} "
        f"({round(len(pcg_df) / len(df) * 100, 2)}%)"
    )
    gene_names = random.choices(list(set(pcg_df["Symbol"])), k=2000)
    nr = 0
    nw = 0
    with GtfIterator("hg38.ncbiRefSeq.gtf") as gtfi, \
            GtfWriter("hg38.ncbiRefSeq_subsampled.gtf") as gtfw:
        for gtf_record in gtfi:
            nr += 1
            gene_id = gtf_record.attribute.get("gene_id")
            if (
                    not gene_id.startswith("TR")
                    and gene_id in gene_names
                    and is_valid_chrname(gtf_record.seqname)
            ):
                nw += 1
                gtfw.write_feature(gtf_record)
    print(
        f"Filtered {nw} gtf records out of {nr} "
        f"({round(nw / nr * 100, 2)}%)"
    )
