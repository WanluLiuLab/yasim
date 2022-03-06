from coverage.annotate import os

import test_tetgs
from bioutils.datastructure.gene_view import GeneView
from commonutils import shell_utils
from commonutils.io.safe_io import get_writer
from commonutils.stdlib_helper import logger_helper

logger_helper.set_level(8)

test_path = test_tetgs.initialize(__name__)

gene_gtf = """
chrI	ncbiRefSeq	exon	4221	4358	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "1"; exon_id "NM_058260.4.1"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	5195	5296	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "2"; exon_id "NM_058260.4.2"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	6037	6327	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "3"; exon_id "NM_058260.4.3"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	9727	9846	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "4"; exon_id "NM_058260.4.4"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	10095	10148	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "5"; exon_id "NM_058260.4.5"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	11641	11689	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "1"; exon_id "NM_058259.4.1"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	14951	15160	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "2"; exon_id "NM_058259.4.2"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	16473	16585	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "3"; exon_id "NM_058259.4.3"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	15103	15160	.	+	.	gene_id "nlp-40"; transcript_id "NM_001306277.1"; exon_number "1"; exon_id "NM_001306277.1.1"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	16473	16585	.	+	.	gene_id "nlp-40"; transcript_id "NM_001306277.1"; exon_number "2"; exon_id "NM_001306277.1.2"; gene_name "nlp-40";
"""


def test_gene() -> None:
    global gene_gtf
    fh = get_writer(os.path.join(test_path, "1.gtf.gz"))
    fh.write(gene_gtf)
    fh.close()
    gv = GeneView.from_file(os.path.join(test_path, "1.gtf.gz"))
    assert list(gv.genes.keys()) == ['homt-1', 'nlp-40']
    assert list(gv.transcripts.keys()) == ['NM_058260.4', 'NM_058259.4', 'NM_001306277.1']
    assert gv.transcripts['NM_058260.4'].exons[0].start == 4221
    gv.to_file(os.path.join(test_path, "2.gtf"))
    gv.standardize()
    gv.to_file(os.path.join(test_path, "3.gtf"))
    os.system(f"gedit {test_path}/*.gtf")
    shell_utils.rm_rf(test_path)
