import json

from yasim_scripts._main.compare_gtts import *


def jsond(dst_json_path: str, tes):
    with open(dst_json_path, "w") as w:
        json.dump(
            {k:list(v) for k, v in tes.items()},
            w
        )

if __name__ == "__main__":
    jsond(
        "aln/ce11.rmsk_consensus.blat.blast9.tes.json",
        convert_aln_blast6_to_tes("aln/ce11.rmsk_consensus.blat.blast9", False),
    )
    jsond(
        "aln/ce11.rmsk_loci.blat.blast9.tes.json",
        convert_aln_blast6_to_tes("aln/ce11.rmsk_loci.blat.blast9", False),
    )
    jsond(
        "aln/rmsk_hmmer.tes.json",
        convert_rmsk_gff_to_tes("aln/rmsk_hmmer.d/ce11_denovo_test.fa.out.gff"),
    )
    jsond(
        "aln/rmsk_rmblast.tes.json",
        convert_rmsk_gff_to_tes("aln/rmsk_rmblast.d/ce11_denovo_test.fa.out.gff"),
    )
    jsond(
        "aln/ce11.minimap2.fc.tes.json",
        convert_fc_assignment_to_tes("aln/ce11.minimap2.fc.tsv.d/ce11.minimap2.bam.featureCounts"),
    )
    jsond(
        "sim/ce11_denovo_test.tes.json",
        convert_translation_instruction_to_tes("sim/ce11_denovo_test.json"),
    )


