import argparse
import random
from typing import List

import pandas as pd

from labw_utils.bioutils.parser.feature import GtfIterator, GtfWriter
from labw_utils.bioutils.algorithm.sequence import is_valid_chrname
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

_lh = get_logger(__name__)


def _parse_args(args: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i',
        '--ncbi_dataset',
        required=True,
        help="The `ncbi_dataset.tsv` retrived from <https://www.ncbi.nlm.nih.gov/labs/data-hub/gene/table/taxon/9606/>",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-g',
        '--gtf',
        required=True,
        help="GTF to be sampled from",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '-o',
        '--out',
        required=True,
        help="Output GTF",
        nargs='?',
        type=str,
        action='store'
    )
    parser.add_argument(
        '--num_genes_to_sample',
        required=False,
        help="Number of genes to be sampled",
        nargs='?',
        type=int,
        default=2000,
        action='store'
    )
    parser.add_argument(
        '--gene_name_attribute',
        required=False,
        help="Attribute name that shows Entrez symbol. Should be `gene_id` or `gene_name`.",
        nargs='?',
        type=str,
        default="gene_id",
        action='store'
    )
    return parser.parse_args(args)


def main(args: List[str]):
    args = _parse_args(args)
    df = pd.read_table(args.ncbi_dataset)
    pcg_df = (
        df.
        query("`Gene Type` == 'PROTEIN_CODING'").
        query("`Gene Group Method` == 'NCBI Ortholog'")
    )
    _lh.info(
        "Filtered %d genes out of %d (%.2f%%)",
        len(pcg_df), len(df), len(pcg_df) / len(df) * 100
    )
    gene_names = random.choices(list(set(pcg_df["Symbol"])), k=args.num_genes_to_sample)
    nr = 0
    nw = 0
    with GtfIterator(args.gtf) as gtfi, \
            GtfWriter(args.out) as gtfw:
        for gtf_record in gtfi:
            nr += 1
            gene_id = gtf_record.attribute.get(args.gene_name_attribute)
            if (
                    not gene_id.startswith("TR")
                    and gene_id in gene_names
                    and is_valid_chrname(gtf_record.seqname)
            ):
                nw += 1
                gtfw.write_feature(gtf_record)
    _lh.info(
        "Filtered %d gtf records out of %d (%.2f%%)",
        nw, nr, nw / nr * 100
    )