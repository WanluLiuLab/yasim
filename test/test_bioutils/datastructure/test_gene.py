import test_tetgs
from bioutils.datastructure.gene_view import GeneView
from commonutils import logger, shell_utils
from commonutils.io.safe_io import get_writer

logger.set_level(8)

test_path = test_tetgs.initialize(__name__)

gene_gtf = """chr1	refGene	transcript	173477335	173488815	.	+	.	gene_id "PRDX6"; transcript_id "NM_004905";  gene_name "PRDX6";
chr1	refGene	exon	173477335	173477492	.	+	.	gene_id "PRDX6"; transcript_id "NM_004905"; exon_number "1"; exon_id "NM_004905.1"; gene_name "PRDX6";
chr1	refGene	5UTR	173477335	173477397	.	+	.	gene_id "PRDX6"; transcript_id "NM_004905"; exon_number "1"; exon_id "NM_004905.1"; gene_name "PRDX6";
chr1	refGene	CDS	173477398	173477492	.	+	0	gene_id "PRDX6"; transcript_id "NM_004905"; exon_number "1"; exon_id "NM_004905.1"; gene_name "PRDX6";
chr1	refGene	exon	173481326	173481482	.	+	.	gene_id "PRDX6"; transcript_id "NM_004905"; exon_number "2"; exon_id "NM_004905.2"; gene_name "PRDX6";
chr1	refGene	CDS	173481326	173481482	.	+	1	gene_id "PRDX6"; transcript_id "NM_004905"; exon_number "2"; exon_id "NM_004905.2"; gene_name "PRDX6";
chr1	refGene	exon	173485361	173485507	.	+	.	gene_id "PRDX6"; transcript_id "NM_004905"; exon_number "3"; exon_id "NM_004905.3"; gene_name "PRDX6";
chr1	refGene	transcript	173477335	173488815	.	+	.	gene_id "PRDX6"; transcript_id "NM_004906";  gene_name "PRDX6";
chr1	refGene	exon	173477335	173477492	.	+	.	gene_id "PRDX6"; transcript_id "NM_004906"; exon_number "1"; exon_id "NM_004906.1"; gene_name "PRDX6";
chr1	refGene	5UTR	173477335	173477397	.	+	.	gene_id "PRDX6"; transcript_id "NM_004906"; exon_number "1"; exon_id "NM_004906.1"; gene_name "PRDX6";
chr1	refGene	CDS	173477398	173477492	.	+	0	gene_id "PRDX6"; transcript_id "NM_004906"; exon_number "1"; exon_id "NM_004906.1"; gene_name "PRDX6";
chr1	refGene	exon	173481326	173481482	.	+	.	gene_id "PRDX6"; transcript_id "NM_004906"; exon_number "2"; exon_id "NM_004906.2"; gene_name "PRDX6";
chr1	refGene	CDS	173481326	173481482	.	+	1	gene_id "PRDX6"; transcript_id "NM_004906"; exon_number "2"; exon_id "NM_004906.2"; gene_name "PRDX6";
chr1	refGene	exon	173485361	173485507	.	+	.	gene_id "PRDX6"; transcript_id "NM_004906"; exon_number "3"; exon_id "NM_004906.3"; gene_name "PRDX6";
chr1	refGene	transcript	173477335	173488815	.	+	.	gene_id "PRDX7"; transcript_id "NM_004907";  gene_name "PRDX7";
chr1	refGene	exon	173477335	173477492	.	+	.	gene_id "PRDX7"; transcript_id "NM_004907"; exon_number "1"; exon_id "NM_004907.1"; gene_name "PRDX7";
chr1	refGene	5UTR	173477335	173477397	.	+	.	gene_id "PRDX7"; transcript_id "NM_004907"; exon_number "1"; exon_id "NM_004907.1"; gene_name "PRDX7";
chr1	refGene	CDS	173477398	173477492	.	+	0	gene_id "PRDX7"; transcript_id "NM_004907"; exon_number "1"; exon_id "NM_004907.1"; gene_name "PRDX7";
chr1	refGene	exon	173481326	173481482	.	+	.	gene_id "PRDX7"; transcript_id "NM_004907"; exon_number "2"; exon_id "NM_004907.2"; gene_name "PRDX7";
chr1	refGene	CDS	173481326	173481482	.	+	1	gene_id "PRDX7"; transcript_id "NM_004907"; exon_number "2"; exon_id "NM_004907.2"; gene_name "PRDX7";
chr1	refGene	exon	173485361	173485507	.	+	.	gene_id "PRDX7"; transcript_id "NM_004907"; exon_number "3"; exon_id "NM_004907.3"; gene_name "PRDX7";
chr1	refGene	transcript	173477335	173488815	.	+	.	gene_id "PRDX7"; transcript_id "NM_004908";  gene_name "PRDX7";
chr1	refGene	exon	173477335	173477492	.	+	.	gene_id "PRDX7"; transcript_id "NM_004908"; exon_number "1"; exon_id "NM_004908.1"; gene_name "PRDX7";
chr1	refGene	5UTR	173477335	173477397	.	+	.	gene_id "PRDX7"; transcript_id "NM_004908"; exon_number "1"; exon_id "NM_004908.1"; gene_name "PRDX7";
chr1	refGene	CDS	173477398	173477492	.	+	0	gene_id "PRDX7"; transcript_id "NM_004908"; exon_number "1"; exon_id "NM_004908.1"; gene_name "PRDX7";
chr1	refGene	exon	173481326	173481482	.	+	.	gene_id "PRDX7"; transcript_id "NM_004908"; exon_number "2"; exon_id "NM_004908.2"; gene_name "PRDX7";
chr1	refGene	CDS	173481326	173481482	.	+	1	gene_id "PRDX7"; transcript_id "NM_004908"; exon_number "2"; exon_id "NM_004908.2"; gene_name "PRDX7";
chr1	refGene	exon	173485361	173485507	.	+	.	gene_id "PRDX7"; transcript_id "NM_004908"; exon_number "3"; exon_id "NM_004908.3"; gene_name "PRDX7";
"""


def test_gene() -> None:
    global gene_gtf
    fh = get_writer(f"{test_path}/1.gtf.gz")
    fh.write(gene_gtf)
    fh.close()
    gv = GeneView.from_file(f"{test_path}/1.gtf.gz")
    assert list(gv.genes.keys()) == ['PRDX6', 'PRDX7']
    assert list(gv.transcripts.keys()) == ['NM_004905', 'NM_004906', 'NM_004907', 'NM_004908']
    assert gv.transcripts['NM_004905'].exons[0].start == 173477335
    gv.to_gtf(f"{test_path}/2.gtf")
    # FIXME: (f"gedit {test_path}/2.gtf")
    shell_utils.rm_rf(test_path)
