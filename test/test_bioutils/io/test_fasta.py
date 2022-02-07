# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: test_fasta.py -- Unit test of corresponding module.
#
#  VERSION HISTORY:
#  2021-08-11 0.1  : Purposed and added by YU Zhejian.
#
# ==============================================================================

"""
test_fasta.py -- Unit test of corresponding module.
"""

import random
import typing

import pytest

import test_tetgs
from bioutils.io import fasta
from commonutils import ioctl, logger

logger.set_level(8)

test_path = test_tetgs.initialize(__name__)

fasta_seq = """>chr1 some att
NNNNNNNNNNNNNNNATCGTTACGTAC
CATATACTATATCTTAGTCTAGTCTAA
CGTCTTTTTCTNNNNNNNNNNNNNNNA
NNNNNNNNATCGTTACGTACTTCTNNN
CATATACTATATCTTAGTCTAGTCTAA
CGTCTTTTTCTNNNNNNNN
>chr2
NNNNNNNNNNNNNNNATCGTTACGTAC
CATATACTATATCTTAGTCTAGTCTAA
CGTCTTTTTCTNNNNNNNNN
>chr3
ACT
ANN
TGN
ATN
ATG
N
"""


def cleanup() -> None:
    ioctl.rm_rf(f"{test_path}/1.fasta.gz")
    ioctl.rm_rf(f"{test_path}/1.fasta")
    ioctl.rm_rf(f"{test_path}/1.fasta.fxi")
    ioctl.rm_rf(f"{test_path}/1.fasta.gz.fai")
    ioctl.rm_rf(f"{test_path}/1.fasta.fai")


def test_rev_compl() -> None:
    assert fasta.reverse_complement("CTGACTGA") == 'TCAGTCAG'


def fasta_with_full_header_assets(fa: fasta.FastaView) -> None:
    assert len(fa) == 3
    assert fa.sequence('chr1 some att', 0, 1) == 'N'
    assert fa.sequence('chr1 some att', 26, 29) == 'CCA'  # Cross the line
    assert fa.sequence('chr1 some att', 28, 29) == 'A'  # Next line
    assert fa.sequence('chr1 some att', 5, 29) == 'NNNNNNNNNNATCGTTACGTACCA'
    with pytest.raises(ValueError):
        fa.sequence('chr1 some att', 5, 1222)
    with pytest.raises(ValueError):
        fa.sequence('chr1 some att', -5, 29)
    with pytest.raises(ValueError):
        fa.sequence('chr1 some att', 500, 29)
    with pytest.raises(ValueError):
        fa.sequence('chr1', 5, 29)
    assert fa.sequence('chr1 some att', 5, 63) == 'NNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTT'
    assert fa.sequence('chr3', 2, 15) == 'TANNTGNATNATG'  # At line end
    assert fa.sequence('chr3', 2, 16) == 'TANNTGNATNATGN'  # Cross the line
    assert fa.sequence('chr2', 0,
                       -1) == 'NNNNNNNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTTCTNNNNNNNNN'


def fasta_without_full_header_assets(fa: fasta.FastaView) -> None:
    assert len(fa) == 3
    assert fa.sequence('chr1', 0, 1) == 'N'
    assert fa.sequence('chr1', 26, 29) == 'CCA'  # Cross the line
    assert fa.sequence('chr1', 28, 29) == 'A'  # Next line
    assert fa.sequence('chr1', 5, 29) == 'NNNNNNNNNNATCGTTACGTACCA'
    with pytest.raises(ValueError):
        fa.sequence('chr1', 5, 1222)
    with pytest.raises(ValueError):
        fa.sequence('chr1', -5, 29)
    with pytest.raises(ValueError):
        fa.sequence('chr1', 500, 29)
    assert fa.sequence('chr1', 5, 63) == 'NNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTT'
    assert fa.sequence('chr3', 2, 15) == 'TANNTGNATNATG'  # At line end
    assert fa.sequence('chr3', 2, 16) == 'TANNTGNATNATGN'  # Cross the line
    assert fa.sequence('chr2', 0,
                       -1) == 'NNNNNNNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTTCTNNNNNNNNN'


def test_fasta_class_without_fai_in_mem() -> None:
    global fasta_seq
    ioctl.rm_rf(f"{test_path}/1.fasta.fai")
    ioctl.rm_rf(f"{test_path}/1.fasta.gz.fai")
    fh = ioctl.get_writer(f"{test_path}/1.fasta.gz")
    fh.write(fasta_seq)
    fh.close()
    fa = fasta.FastaView(f"{test_path}/1.fasta.gz", read_into_memory=True)
    fasta_without_full_header_assets(fa)
    fa.close()
    fa = fasta.FastaView(f"{test_path}/1.fasta.gz", read_into_memory=True, all_header=True)
    fasta_with_full_header_assets(fa)
    fa.close()
    cleanup()


def test_fasta_class_without_fai_without_mem() -> None:
    global fasta_seq
    ioctl.rm_rf(f"{test_path}/1.fasta.fai")
    ioctl.rm_rf(f"{test_path}/1.fasta.gz.fai")
    fh = ioctl.get_writer(f"{test_path}/1.fasta")
    fh.write(fasta_seq)
    fh.close()
    fa = fasta.FastaView(f"{test_path}/1.fasta", read_into_memory=False)
    fasta_without_full_header_assets(fa)
    fa.close()
    fa = fasta.FastaView(f"{test_path}/1.fasta", read_into_memory=False, all_header=True)
    fasta_with_full_header_assets(fa)
    fa.close()
    cleanup()


# noinspection all
def test_fasta_class_with_fai_without_mem() -> None:
    try:
        import pysam
    except ImportError:
        return
    global fasta_seq
    fh = ioctl.get_writer(f"{test_path}/1.fasta")
    fh.write(fasta_seq)
    fh.close()
    pysam.faidx(f"{test_path}/1.fasta")
    fa = fasta.FastaView(f"{test_path}/1.fasta", read_into_memory=False)
    fasta_without_full_header_assets(fa)
    fa.close()
    fa = fasta.FastaView(f"{test_path}/1.fasta", read_into_memory=False, all_header=True)
    with pytest.raises(ValueError):
        fasta_with_full_header_assets(fa)
    fa.close()
    cleanup()


def _fatsa_gen(writer: typing.TextIO) -> []:
    retv = []
    this_tras = str.maketrans('01234', 'NAGCT')

    for i in range(0, 20):
        writer.write(f'>chr{i}')
        chr_len = random.randint(100, 100000)
        chr_len_b = random.randint(10, 100)
        in_seq = ''
        for _ in range(chr_len):
            in_seq = in_seq + str(random.randint(0, 4))
        rp = 0
        while rp < len(in_seq):
            in_seq = in_seq[:rp] + '\n' + in_seq[rp:]
            rp += chr_len_b
        writer.write(in_seq.translate(this_tras) + '\n')
        # chr_len, start, end
        start = 10
        end = 0
        while start > end:
            start = random.randint(0, chr_len)
            end = random.randint(0, chr_len)
        retv.append([chr_len, start, end])
    return retv


def test_dynamic_asserts() -> None:
    try:
        import pyfaidx
    except ImportError:
        return
    fh = ioctl.get_writer(f"{test_path}/2.fasta")
    retv = _fatsa_gen(fh)
    fh.close()
    tru = pyfaidx.FastaView(f"{test_path}/2.fasta", one_based_attributes=False)
    this = fasta.FastaView(f"{test_path}/2.fasta", read_into_memory=False)
    i = 0
    for item in retv:
        assert tru[f'chr{i}'][item[1]: item[2]].seq == this.sequence(f'chr{i}', item[1], item[2])
        i += 1
    tru.close()
    this.close()
    ioctl.rm_rf(f"{test_path}/2.fasta")
    ioctl.rm_rf(f"{test_path}/2.fasta.fai")


if __name__ == "__main__":
    test_rev_compl()
    test_fasta_class_without_fai_in_mem()
    test_fasta_class_without_fai_without_mem()
    test_fasta_class_with_fai_without_mem()
    test_dynamic_asserts()
