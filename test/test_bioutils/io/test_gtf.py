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

import test_tetgs
from commonutils.io.safe_io import get_writer

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

with get_writer(fasta_path) as writer:
    writer.write(fasta_contents)

with get_writer(gtf_path) as writer:
    writer.write(gtf_contents)
