import os

from yasim.helper.tcr import create_tcr_cache



if __name__ == "__main__":

    basename = "."
    create_tcr_cache(
        ref_fa_path=os.path.join(basename, "hg38.fa"),
        ref_gtf_path=os.path.join(basename, "hg38.ncbiRefSeq_chr7_14.gtf"),
        tcr_genelist_path=os.path.join(basename, "tcr_genelist.json"),
        tcr_aa_table_path=os.path.join(basename, "IMGT_Protein_Display.json"),
        tcr_cache_path=os.path.join(basename, "tcr_cache.json")
    )
