import bz2
import gzip
import io
import lzma
import os
import random
import string

import test_tetgs
from commonutils import shell_utils
from commonutils.io import ArchiveBaseIO
from commonutils.io.safe_io import get_writer, get_reader
from commonutils.io.tqdm_reader import get_tqdm_reader, get_tqdm_line_reader

test_path = test_tetgs.initialize(__name__)

available_chars = list(string.ascii_letters + string.digits + "\n")

contents = "".join((random.choice(available_chars) for _ in range(1024 * 1024)))
len_contents = len(contents)

contents_list = contents.splitlines()


def assess_archive_io(filename: str):
    shell_utils.rm_rf(filename)
    with get_writer(filename) as writer:
        writer.write(contents)
    with get_reader(filename) as reader:
        assert reader.read(len_contents) == contents
    with get_tqdm_reader(filename) as reader:
        assert reader.read(len_contents) == contents
    with get_tqdm_line_reader(filename) as reader:
        i = 0
        assert reader._tqdm.total == 1 + len(contents_list)
        for line in reader:
            assert contents_list[i] == line
            i += 1
            assert reader._tqdm._n == i



def test_txt():
    filename = os.path.join(test_path, f"1.txt")
    assess_archive_io(filename)
    bare_archive_io = ArchiveBaseIO(open(filename, "rt"))
    assert bare_archive_io.read(len_contents) == contents
    bare_archive_io.close()
    shell_utils.rm_rf(filename)


def test_bz2():
    filename = os.path.join(test_path, f"1.bz2")
    assess_archive_io(filename)
    bare_archive_io = ArchiveBaseIO(bz2.open(filename, "rt"))
    assert bare_archive_io.read(len_contents) == contents
    bare_archive_io.close()
    shell_utils.rm_rf(filename)


def test_gz():
    filename = os.path.join(test_path, f"1.gz")
    assess_archive_io(filename)
    bare_archive_io = ArchiveBaseIO(gzip.open(filename, "rt"))
    assert bare_archive_io.read(len_contents) == contents
    bare_archive_io.close()
    shell_utils.rm_rf(filename)


def test_xz():
    filename = os.path.join(test_path, f"1.xz")
    assess_archive_io(filename)
    bare_archive_io = ArchiveBaseIO(lzma.open(filename, "rt"))
    assert bare_archive_io.read(len_contents) == contents
    bare_archive_io.close()
    shell_utils.rm_rf(filename)


def test_lzma():
    filename = os.path.join(test_path, f"1.lzma")
    bare_archive_io = ArchiveBaseIO(lzma.open(filename, "rt"))
    assert bare_archive_io.read(len_contents) == contents
    bare_archive_io.close()
    shell_utils.rm_rf(filename)

def test_string_io():
    # FIXME
    sio = io.StringIO(contents)
    bare_archive_io = ArchiveBaseIO(sio)
    assert bare_archive_io.read(len_contents) == contents


