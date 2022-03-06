"""
test_fasta.py -- Unit test of corresponding module.
"""
import os
import random
from typing import List, IO, Tuple

import pytest

import test_tetgs
from bioutils.datastructure.fasta_view import FastaView, create_fai
from commonutils import shell_utils
from commonutils.io.safe_io import get_writer
from commonutils.stdlib_helper import logger_helper

logger_helper.set_level(8)

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
>chr4
AAAAAAAAAACCCCCC
>chr6
CTA
"""

VALID_NEWLINE = ("\n", "\r\n")

FASTA_FILENAME = os.path.join(test_path, "1.fasta")
FASTA_FAI_FILENAME = FASTA_FILENAME + ".fai"


def cleanup() -> None:
    shell_utils.rm_rf(FASTA_FAI_FILENAME)
    shell_utils.rm_rf(FASTA_FILENAME)


def _public_asserts(fa: FastaView) -> None:
    print(fa.chr_names)
    assert len(fa) == 5
    assert fa.sequence('chr3', 2, 15) == 'TANNTGNATNATG'  # At line end
    assert fa.sequence('chr3', 2, 16) == 'TANNTGNATNATGN'  # Cross the line
    assert fa.sequence('chr2') == 'NNNNNNNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTTCTNNNNNNNNN'
    assert fa.sequence('chr2') == fa.sequence('chr2', 0)
    assert fa.sequence('chr2') == fa.sequence('chr2', 0, -1)
    assert fa.sequence('chr4') == 'AAAAAAAAAACCCCCC'
    assert fa.sequence('chr6') == 'CTA'
    with pytest.raises(ValueError):
        fa.sequence('chr2', 5, 1222)
    with pytest.raises(ValueError):
        fa.sequence('chr2', -5, 29)
    with pytest.raises(ValueError):
        fa.sequence('chr2', 500, 29)


def _fasta_with_full_header_assets(fa: FastaView) -> None:
    assert fa.sequence('chr1 some att', 0, 1) == 'N'
    assert fa.sequence('chr1 some att', 26, 29) == 'CCA'  # Cross the line
    assert fa.sequence('chr1 some att', 28, 29) == 'A'  # Next line
    assert fa.sequence('chr1 some att', 5, 29) == 'NNNNNNNNNNATCGTTACGTACCA'
    with pytest.raises(ValueError):
        fa.sequence('chr1', 5, 29)
    assert fa.sequence('chr1 some att', 5, 63) == 'NNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTT'
    _public_asserts(fa)


def _fasta_without_full_header_assets(fa: FastaView) -> None:
    assert fa.sequence('chr1', 0, 1) == 'N'
    assert fa.sequence('chr1', 26, 29) == 'CCA'  # Cross the line
    assert fa.sequence('chr1', 28, 29) == 'A'  # Next line
    assert fa.sequence('chr1', 5, 29) == 'NNNNNNNNNNATCGTTACGTACCA'
    with pytest.raises(ValueError):
        fa.sequence('chr1 some att', 5, 29)
    assert fa.sequence('chr1', 5, 63) == 'NNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTT'
    _public_asserts(fa)


def test_fasta_class_without_fai_in_mem() -> None:
    global fasta_seq
    cleanup()
    for newline in VALID_NEWLINE:
        fh = get_writer(FASTA_FILENAME, newline=newline)
        fh.write(fasta_seq)
        fh.close()
        fa = FastaView(FASTA_FILENAME, read_into_memory=True, full_header=False)
        _fasta_without_full_header_assets(fa)
        fa.close()
        fa = FastaView(FASTA_FILENAME, read_into_memory=True, full_header=True)
        _fasta_with_full_header_assets(fa)
        fa.close()
        cleanup()


def test_fasta_class_without_fai_without_mem() -> None:
    global fasta_seq
    cleanup()
    for newline in VALID_NEWLINE:
        fh = get_writer(FASTA_FILENAME, newline=newline)
        fh.write(fasta_seq)
        fh.close()
        fa = FastaView(FASTA_FILENAME, read_into_memory=False)
        _fasta_without_full_header_assets(fa)
        fa.close()
        shell_utils.rm_rf(FASTA_FAI_FILENAME)
        fa = FastaView(FASTA_FILENAME, read_into_memory=False, full_header=True)
        with pytest.raises(ValueError):
            _fasta_with_full_header_assets(fa)
        fa.close()
        cleanup()


def test_fai() -> None:
    try:
        from pysam import faidx
    except ImportError:
        return
    global fasta_seq
    for newline in VALID_NEWLINE:
        fh = get_writer(FASTA_FILENAME, newline=newline)
        fh.write(fasta_seq)
        fh.close()
        faidx(FASTA_FILENAME)
        create_fai(FASTA_FILENAME, os.path.join(test_path, "1_tetgs.fai"))
        assert open(FASTA_FAI_FILENAME).read() == open(os.path.join(test_path, "1_tetgs.fai")).read()
        cleanup()


def test_fasta_class_with_fai_without_mem() -> None:
    try:
        from pysam import faidx
    except ImportError:
        return
    global fasta_seq
    for newline in VALID_NEWLINE:
        print(len(newline))
        fh = get_writer(FASTA_FILENAME, newline=newline)
        fh.write(fasta_seq)
        fh.close()
        faidx(FASTA_FILENAME)
        fa = FastaView(FASTA_FILENAME, read_into_memory=False)
        _fasta_without_full_header_assets(fa)
        fa.close()
        fa = FastaView(FASTA_FILENAME, read_into_memory=False, full_header=True)
        with pytest.raises(ValueError):
            _fasta_with_full_header_assets(fa)
        fa.close()
        cleanup()


def _fatsa_gen(writer: IO) -> List[Tuple[int, int, int]]:
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
        retv.append((chr_len, start, end))
    return retv


def _dynamic_asserts(retv):
    try:
        import pyfaidx
    except ImportError:
        return
    tru = pyfaidx.FastaView(f"{test_path}/2.fasta", one_based_attributes=False)
    for this in (
            FastaView(f"{test_path}/2.fasta", read_into_memory=False),
            FastaView(f"{test_path}/2.fasta", read_into_memory=True)
    ):
        i = 0
        for item in retv:
            assert tru[f'chr{i}'][item[1]: item[2]].seq == this.sequence(f'chr{i}', item[1], item[2])
            i += 1
        this.close()
    tru.close()
    shell_utils.rm_rf(f"{test_path}/2.fasta")
    shell_utils.rm_rf(f"{test_path}/2.fasta.fai")


def test_dynamic_asserts() -> None:
    try:
        import pyfaidx
    except ImportError:
        return
    for newline in ('\n', '\r\n'):
        fh = get_writer(f"{test_path}/2.fasta", newline=newline)
        retv = _fatsa_gen(fh)
        fh.close()
        _dynamic_asserts(retv)


if __name__ == "__main__":
    test_fasta_class_without_fai_in_mem()
    test_fasta_class_without_fai_without_mem()
    test_fasta_class_with_fai_without_mem()
    test_dynamic_asserts()
