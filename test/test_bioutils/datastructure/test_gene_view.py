import os

import pytest

import conftest
from bioutils.datastructure.gene_view import GeneViewFactory
from commonutils import shell_utils
from commonutils.io.safe_io import get_writer

gene_gtf = """
chrI	ncbiRefSeq	exon	4221	4358	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "1"; exon_id "NM_058260.4.1"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	5195	5296	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "2"; exon_id "NM_058260.4.2"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	6037	6327	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "3"; exon_id "NM_058260.4.3"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	9727	9846	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "4"; exon_id "NM_058260.4.4"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	10095	10148	.	-	.	gene_id "homt-1"; transcript_id "NM_058260.4"; exon_number "5"; exon_id "NM_058260.4.5"; gene_name "homt-1";
chrI	ncbiRefSeq	exon	11641	11689	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "1"; exon_id "NM_058259.4.1"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	11641	11689	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "1"; exon_id "NM_058259.4.1"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	16473	16585	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "3"; exon_id "NM_058259.4.3"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	14951	15160	.	+	.	gene_id "nlp-40"; transcript_id "NM_058259.4"; exon_number "2"; exon_id "NM_058259.4.2"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	15103	15160	.	+	.	gene_id "nlp-40"; transcript_id "NM_001306277.1"; exon_number "1"; exon_id "NM_001306277.1.1"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	16473	16585	.	+	.	gene_id "nlp-40"; transcript_id "NM_001306277.1"; exon_number "2"; exon_id "NM_001306277.1.2"; gene_name "nlp-40";
chrI	ncbiRefSeq	exon	8484719	8484913	.	-	.	gene_id "D1081.6"; transcript_id "NM_059899.3"; exon_number "1"; exon_id "NM_059899.3.1"; gene_name "D1081.6";
chrI	ncbiRefSeq	transcript	6492885	6493553	0	-	.	gene_id "mdt-18"; transcript_id "NM_001322685.1"; gene_name "mdt-18"; 
"""


@pytest.fixture(scope="module", autouse=True)
def initialize_module(initialize_session) -> conftest.ModuleTestInfo:
    """
    This function sets up a directory for testing
    """
    session_test_info = initialize_session
    module_test_info = conftest.ModuleTestInfo(session_test_info.base_test_dir, __name__)
    with get_writer(os.path.join(module_test_info.path, "1.gtf")) as fh:
        fh.write(gene_gtf)
    yield module_test_info
    module_test_info.teardown()


def test_gene(initialize_module) -> None:
    test_path = initialize_module.path
    gv = GeneViewFactory.from_file(os.path.join(test_path, "1.gtf"))
    assert list(gv.genes.keys()) == ['homt-1', 'nlp-40', 'D1081.6', "mdt-18"]
    assert list(gv.transcripts.keys()) == ['NM_058260.4', 'NM_058259.4', 'NM_001306277.1', 'NM_059899.3',
                                           "NM_001322685.1"]
    assert gv.transcripts['NM_058260.4'].exons[0].start == 4221
    gv.standardize()
    assert list(gv.transcripts.keys()) == ['NM_058260.4', 'NM_058259.4', 'NM_001306277.1', 'NM_059899.3']
    gv.to_file(os.path.join(test_path, "3.gtf"))
    os.system(f"gedit {test_path}/*.gtf")
    shell_utils.rm_rf(test_path)
