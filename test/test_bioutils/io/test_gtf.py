# ==============================================================================
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: test_gtf.py -- Unit test of corresponding module.
#
#  VERSION HISTORY:
#  2021-08-16 0.1  : Added by YU Zhejian.
#
# ==============================================================================
"""
test_gtf.py -- Unit test of corresponding module.
"""
import time

import test_tetgs
from bioutils.io import gtf
from bioutils.io.gtf import GtfRecord, parse_gtf_attrs
from commonutils import ioctl
from commonutils.str_utils import to_dict

test_path = test_tetgs.initialize(__name__)

# The following file was retrieved from <http://genome.ucsc.edu/cgi-bin/hgTables>
gtf_file_usuc = """# geno_name	?	?	geno_start	geno_end	swScore	strand	?	rep_name
# No idea how the transcript ID is generated.
chr1	hg38_rmsk	exon	67108754	67109046	1892.000000	+	.	gene_id "L1P5"; transcript_id "L1P5"; 
chr1	hg38_rmsk	exon	8388316	8388618	2582.000000	-	.	gene_id "AluY"; transcript_id "AluY"; 
chr1	hg38_rmsk	exon	25165804	25166380	4085.000000	+	.	gene_id "L1MB5"; transcript_id "L1MB5"; 
chr1	hg38_rmsk	exon	33554186	33554483	2285.000000	-	.	gene_id "AluSc"; transcript_id "AluSc"; 
chr1	hg38_rmsk	exon	41942895	41943205	2451.000000	-	.	gene_id "AluY"; transcript_id "AluY_dup1"; 
chr1	hg38_rmsk	exon	50331337	50332274	1587.000000	+	.	gene_id "HAL1"; transcript_id "HAL1"; 
"""

fasta_contents = """>chr1
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACAGTACTGGCGGATTATAGGGAAACACCCGGC
AGCATATGCTGTTTGGTCTTTGCTGTTCCTGCATGTAGTTTAAACGAGATTGCCAGCACCGGGTATCATTCACCATTTTC
TCTTTTCGTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTTCTGTTCATGTGTATTTGCTGTCTCTTAGCCCAGAC
CTTCCCGTGTCCTTTCCACCGGGCCTTTGAGAGGTCACAGGGTCTTGATGCTGTGGTCTTCATCTGCAGGTGTCTGACTC
TCCAGCAACTGCTGGCCTGTGCCAGGGTGCAAGCTGAGCACTGGAGTGGAGTTTTCCTGTGGAGAGGAGCCATGCCTAGC
AGTGGGATGGGCCATTGTTCATCTTCTGGCCCCTGTTGTCTGCATGTAACTTAATACCACAACCAGGCATAGGGGAAAGC
ggccaggagcagtggctcacgcctgtaatcccagcactttaggaggccggggcgggcagatcacaaagtcaggagatcgc
agaccatcctagccaacatggtgaaaccccgtctctactaaaaacacaaaaattagctgggcatggtggcacgcgcctgc
taatcccagctacttatgaggctgaggcaggagaatcgcttgaacccaggaggcggaggttgcagtgagccgagatcgcc
accactgccctccagcctggcgacagagtgagactccatctcaaaaaaaaataaaaaATTGGAGGAAAGATGAGTGAGAC
GCATCAACTTCTCTCACAACCTAGGCCAGTAAGTAGTGCTTGTGCTCATCTCCTTGGCTGTGATACGTGGCCGGCCCTCC
GCTCCAGCAGCTGGACCCCTACCTGCCGTCTGCTGCCATCGGAGCCCAAAGCCGGGCTGTGACTGCTCAGACCAGCCGGC
CTGGAGGGAGGGGCTCAGCAGGTCTGGCTTTGGCCCTGGGAGAGCAGGTGGAAGATCAGGCAGGCCATCGCTGCCACAGC
AACCCAGTGGATTGGCCTAGGTGGGATCTCTGAGCTCAACAAGCCCTCTCTGGGTGGTAGGTGCAGAGACGGGAGGGGCC
AGAGCCGCAGGCACAGCCAAGAGGGCTGAAGAAATGGTAGAACGGAGCAGCTGGTGATGTGTGGGCCCACCGGCCCCAGC
GCTCCTGTCTCCCCCCAGGTGTGTGGTGATGCCAGGCATGCCCTTCCCCAGCATCAGGTCTCCAGAGCTGCAGAAGACGC
ACGGCCGACTTGGATCACACTCTTGTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
"""

gtf_contents = """
chr1	hg38_rmsk	exon	481	777	1892.000000	+	.	gene_id "MYTE1"; transcript_id "MYTE1_dup0"; 
"""

fasta_path = f'{test_path}/2.fa'
gtf_path = f'{test_path}/2.gtf'

with ioctl.get_writer(fasta_path) as writer:
    writer.write(fasta_contents)

with ioctl.get_writer(gtf_path) as writer:
    writer.write(gtf_contents)

def test_time():
    gtf_str = 'chr1	hg38_rmsk	exon	50331337	50332274	1587.000000	+	.	gene_id "HAL1"; transcript_id "HAL1"; '
    try:
        from bioutils.io.gtf._gtf_record import GtfRecord as pyGtfRecord
        from bioutils.io.gtf._c_gtf_record import GtfRecord as cGtfRecord
    except ImportError:
        return
    time_start = time.time()
    times = 100000
    for i in range(times):
        _ = pyGtfRecord.from_string(gtf_str)
    time_py = time.time()
    for j in range(times):
        _ = cGtfRecord.from_string(gtf_str)
    time_c = time.time()
    print(f"time_py={time_py-time_start} time_c={time_c - time_py}")


def test_attributes_parser_time():
    time_start = time.time()
    times = 1000000
    for i in range(times):
        _ = parse_gtf_attrs(
        b'gene_id "STRG.23"; transcript_id "STRG.23.1"; reference_id "XM_017000355.2"; ref_gene_id "MIB2"; ref_gene_name "MIB2"; cov "1.073021"; FPKM "196.106232"; TPM "586.075928";')
    time_parse_gtf_attrs = time.time()
    for j in range(times):
        _ = to_dict(
        'gene_id "STRG.23"; transcript_id "STRG.23.1"; reference_id "XM_017000355.2"; ref_gene_id "MIB2"; ref_gene_name "MIB2"; cov "1.073021"; FPKM "196.106232"; TPM "586.075928";',
            field_sep = ' ', record_sep=';', resolve_str=True, quotation_mark='\'\"'
        )
    time_to_dict = time.time()
    print(f"time_parse_gtf_attrs={time_parse_gtf_attrs-time_start} time_to_dict={time_to_dict - time_parse_gtf_attrs}")


def test_base():
    gtf_obj = gtf.SimpleGtfView(gtf_path)
    gtf_obj.attach_external_fasta(fasta_path)
    a = list(gtf_obj.fetch('chr1'))[0]

    assert parse_gtf_attrs(
        b'gene_id "STRG.23"; transcript_id "STRG.23.1"; reference_id "XM_017000355.2"; ref_gene_id "MIB2"; ref_gene_name "MIB2"; cov "1.073021"; FPKM "196.106232"; TPM "586.075928";') == {'gene_id': 'STRG.23', 'transcript_id': 'STRG.23.1', 'reference_id': 'XM_017000355.2', 'ref_gene_id': 'MIB2', 'ref_gene_name': 'MIB2', 'cov': 1.073021, 'FPKM': 196.106232}

    assert a.sequence == 'ggccaggagcagtggctcacgcctgtaatcccagcactttaggaggccggggcgggcagatcacaaagtcaggagatcgcagaccatcctagccaacatggtgaaaccccgtctctactaaaaacacaaaaattagctgggcatggtggcacgcgcctgctaatcccagctacttatgaggctgaggcaggagaatcgcttgaacccaggaggcggaggttgcagtgagccgagatcgccaccactgccctccagcctggcgacagagtgagactccatctcaaaaaaaaataaaaa'


# def test_fetch_time():
#     # FIXME: File do not exist
#     # FIXME: Need time specifications
#     gtf_obj = gtf.SimpleGtfView('/home/yuzj/Desktop/BioRef/GRCh38_GENCODE_rmsk_TE_chr1.gtf')
#     for _ in range(1):
#         start = random.randint(1, 503322730)
#         NANOPORE_MAX_READ_LENGTH = 2000000
#         end = random.randint(start, min(503322740, start + NANOPORE_MAX_READ_LENGTH))
#         _ = gtf_obj.fetch('chr1', start, end)
#         assert _ is not None


def test_gtf_record():
    gtf_str = 'chr1	hg38_rmsk	exon	50331337	50332274	1587.000000	+	.	gene_id "HAL1"; transcript_id "HAL1"; '
    gtf_from_line = GtfRecord.from_string(gtf_str)
    assert gtf_from_line.seqname == 'chr1'
    assert gtf_from_line.source == 'hg38_rmsk'
    assert gtf_from_line.feature == 'exon'
    assert gtf_from_line.start == 50331337
    assert gtf_from_line.end == 50332274
    assert gtf_from_line.score == 1587
    assert gtf_from_line.strand == '+'
    assert gtf_from_line.attribute['gene_id'] == 'HAL1'
    assert str(
        gtf_from_line) == 'chr1	hg38_rmsk	exon	50331337	50332274	1587.0	+	.	gene_id "HAL1"; transcript_id "HAL1"; '


def test_gtf_reader():
    global gtf_file_usuc
    with ioctl.get_writer(f"{test_path}/1.gtf") as writer:
        writer.write(gtf_file_usuc)
    gtf_file = gtf.SimpleGtfView(f"{test_path}/1.gtf")
    assert len(gtf_file) == 6
    gtf_file.add_record(GtfRecord.from_string(
        'chr1	hg38_rmsk	exon	50331337	50332274	1587	+	.	gene_id "HAL1"; transcript_id "HAL2"; '))
    assert len(gtf_file) == 7
    gtf_file.to_file(f"{test_path}/1.gtf")
    del gtf_file
    gtf_file = gtf.SimpleGtfView(f"{test_path}/1.gtf")
    assert len(gtf_file) == 7
    gtf_file = gtf.SimpleGtfView(f"{test_path}/1.gtf")
    assert len(gtf_file) == 7
    del gtf_file
    ioctl.rm_rf(f"{test_path}/1.gtf")


def test_gtf_with_intervaltree():
    """
    Interval-tree specific tests.
    """
    with ioctl.get_writer(f"{test_path}/1.gtf") as writer:
        writer.write(gtf_file_usuc)
    gtf_file = gtf.SimpleGtfView(f"{test_path}/1.gtf")
    assert list(gtf_file.fetch('chr7', 50331337, 50332274)) == []
    assert len(list(gtf_file.fetch('chr1', 50331337, 50332274))) == 1
    assert len(list(gtf_file.fetch('chr1', 50332273, 50332275))) == 1
    assert len(list(gtf_file.fetch('chr1', 50331336, 50331338))) == 1
    assert len(list(gtf_file.fetch('chr1', 50331333, 50332276))) == 1
    assert len(list(gtf_file.fetch('chr1', 50331344, 50332264))) == 1
    assert len(list(gtf_file.fetch('chr1', 50332275, 50332276))) == 0
    assert len(list(gtf_file.fetch('chr1', 50331335, 50331336))) == 0
    assert len(list(gtf_file.fetch('chr1', 41942915, 50332266))) == 2



if __name__ == "__main__":
    test_gtf_reader()
    test_gtf_record()
    test_gtf_with_intervaltree()
