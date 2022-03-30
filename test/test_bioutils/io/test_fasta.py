"""
test_fasta.py -- Unit test of corresponding module.
"""
import os

import pytest

import conftest
from bioutils.datastructure.fasta_view import FastaViewFactory, FastaViewType
from bioutils.io.fai import create_fai_from_fasta
from commonutils import shell_utils
from commonutils.io.safe_io import get_writer
from commonutils.stdlib_helper import logger_helper

lh = logger_helper.get_logger(__name__)

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


class TestFastaModuleInfo(conftest.ModuleTestInfo):
    fasta_filename: str
    fasta_index_filename: str

    def __init__(self, base_test_dir: str, name: str):
        super().__init__(base_test_dir, name)
        self.fasta_filename = os.path.join(self.path, "1.fasta")
        self.fasta_index_filename = self.path + ".fai"

    def teardown(self):
        self.cleanup_intermediate_files()
        super().teardown()

    def cleanup_intermediate_files(self):
        shell_utils.rm_rf(self.fasta_filename)
        shell_utils.rm_rf(self.fasta_index_filename)


@pytest.fixture(scope="module", autouse=True)
def initialize_module(initialize_session) -> TestFastaModuleInfo:
    """
    This function sets up a directory for testing
    """
    session_test_info = initialize_session
    module_test_info = TestFastaModuleInfo(session_test_info.base_test_dir, __name__)
    yield module_test_info
    module_test_info.teardown()


VALID_NEWLINE = ("\n", "\r\n")


def _public_asserts(fa: FastaViewType) -> None:
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


def _fasta_with_full_header_assets(fa: FastaViewType) -> None:
    assert fa.sequence('chr1 some att', 0, 1) == 'N'
    assert fa.sequence('chr1 some att', 26, 29) == 'CCA'  # Cross the line
    assert fa.sequence('chr1 some att', 28, 29) == 'A'  # Next line
    assert fa.sequence('chr1 some att', 5, 29) == 'NNNNNNNNNNATCGTTACGTACCA'
    with pytest.raises(ValueError):
        fa.sequence('chr1', 5, 29)
    assert fa.sequence('chr1 some att', 5, 63) == 'NNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTT'
    _public_asserts(fa)


def _fasta_without_full_header_assets(fa: FastaViewType) -> None:
    assert fa.sequence('chr1', 0, 1) == 'N'
    assert fa.sequence('chr1', 26, 29) == 'CCA'  # Cross the line
    assert fa.sequence('chr1', 28, 29) == 'A'  # Next line
    assert fa.sequence('chr1', 5, 29) == 'NNNNNNNNNNATCGTTACGTACCA'
    with pytest.raises(ValueError):
        fa.sequence('chr1 some att', 5, 29)
    assert fa.sequence('chr1', 5, 63) == 'NNNNNNNNNNATCGTTACGTACCATATACTATATCTTAGTCTAGTCTAACGTCTTTTT'
    _public_asserts(fa)


@pytest.mark.parametrize(
    argnames="kwargs",
    argvalues=(
            {"read_into_memory": True},
            {"read_into_memory": False}
    ),
    ids=(
            "test_read_into_memory",
            "test_no_read_into_memory"
    )
)
def test_newline_with_or_without_full_header(initialize_module, kwargs) -> None:
    for newline in VALID_NEWLINE:
        with get_writer(initialize_module.fasta_filename, newline=newline) as fh:
            fh.write(fasta_seq)
        with FastaViewFactory(initialize_module.fasta_filename, full_header=False, **kwargs) as fa:
            _fasta_without_full_header_assets(fa)
        if kwargs['read_into_memory']:
            with FastaViewFactory(initialize_module.fasta_filename, full_header=True, **kwargs) as fa:
                _fasta_with_full_header_assets(fa)
        initialize_module.cleanup_intermediate_files()


def test_fai(initialize_module) -> None:
    try:
        from pysam import faidx
    except ImportError:
        return
    global fasta_seq
    for newline in VALID_NEWLINE:
        fasta_filename = initialize_module.fasta_filename
        with get_writer(fasta_filename, newline=newline) as fh:
            fh.write(fasta_seq)
        faidx(fasta_filename)
        create_fai_from_fasta(fasta_filename, fasta_filename + ".tetgs.fai")
        assert open(fasta_filename + ".fai").read() == open(fasta_filename + ".tetgs.fai").read()
        initialize_module.cleanup_intermediate_files()
