import rmsk_parser
from labw_utils.bioutils.parser.gtf import GtfIteratorWriter

RMSK_OUT_GFF_PATH = "/slurm/home/yrd/liulab/yuzhejian/yasim/explore/test_te_quant/ref/rmsk_wd/ce11.out.gff"
OUT_GTF_PATH = "/slurm/home/yrd/liulab/yuzhejian/yasim/explore/test_te_quant/ref/rmsk_wd/ce11.out.filtered.gtf"

with GtfIteratorWriter(OUT_GTF_PATH) as gtfw:
    with rmsk_parser.RMSKGffIterator(RMSK_OUT_GFF_PATH, True) as rmski:
        for i, feature in enumerate(rmski):
            repeat_name = feature.attribute_get("repeat_name")
            if feature.end0b - feature.start0b + 1 < 20 or "(" in repeat_name or repeat_name.endswith("-rich"):
                continue
            gtfw.write(
                feature.update_attribute(gene_id=repeat_name, transcript_id=repeat_name + f"_{i}").update(
                    feature=repeat_name + f"{i}"
                )
            )
