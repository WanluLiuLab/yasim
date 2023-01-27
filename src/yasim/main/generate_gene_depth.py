from typing import List

from labw_utils.bioutils.datastructure.gene_view import GeneViewFactory

from yasim.helper import depth
from yasim.main.generate_depth_v2 import _parse_args


def main(args: List[str]):
    args = _parse_args(args)
    gv = GeneViewFactory.from_file(args.gtf)
    dge_data = depth.simulate_gene_level_depth_gmm(gv=gv, mu=args.mu)
    depth.write_depth(dge_data, args.out, feature_name="GENE_ID")
