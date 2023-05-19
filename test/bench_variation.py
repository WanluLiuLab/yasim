import numpy as np

from labw_utils.mlutils.ndarray_helper import describe
from yasim.helper.depth import simulate_gene_level_depth_gmm, simulate_isoform_variance_inside_a_gene

from labw_utils.bioutils.datastructure.gene_view_v0_1_x.gene_view import GeneViewFactory
from yasim.helper.plot_utils import plot

if __name__ == '__main__':
    gv = GeneViewFactory.from_file("ce11.ncbiRefSeq.gtf")
    for _ in range(20):
        d = simulate_gene_level_depth_gmm(
            gv=gv,
            mu=40,
            low_cutoff=1,
            high_cutoff_ratio=10E4
        )
        gene_expr_levels = np.array(list(d.values()), dtype="float")
        inside_isoform_vars = []
        isoform_expr_levels = []
        for gene_name, mean_expr in d.items():
            i = simulate_isoform_variance_inside_a_gene(
                n=gv.get_gene(gene_name).number_of_transcripts,
                mu=mean_expr,
                alpha=10,
                low_cutoff=1,
                high_cutoff_ratio=10E4
            )
            isoform_expr_levels.extend(i)
            inside_isoform_vars.append(max(i) / min(i))
        isoform_expr_levels = np.array(isoform_expr_levels)
        inside_isoform_vars = np.array(inside_isoform_vars)
        print("gene_expr_levels: " + describe(gene_expr_levels))
        print("inside_isoform_vars: " + describe(inside_isoform_vars))
        print("isoform_expr_levels: " + describe(isoform_expr_levels))
        # plot(np.log(gene_expr_levels), title="log(gene_expr_levels)")
        # plot(np.log(isoform_expr_levels), title="log(isoform_expr_levels)")
