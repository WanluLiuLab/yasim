# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: test_io.py -- Unit test of corresponding module.
#
#  VERSION HISTORY:
#  2021-08-10 0.1  : Added by YU Zhejian.
#
# ==============================================================================
"""
test_io.py -- Unit test of corresponding module.
"""
import random
import string

import test_tetgs
from commonutils import shell_utils
from commonutils.io.safe_io import get_writer, get_reader
from commonutils.io.tqdm_reader import get_tqdm_reader, get_tqdm_line_reader

test_path = test_tetgs.initialize(__name__)

available_chars = list(string.ascii_letters + string.digits + "\n")

contents = "".join((random.choice(available_chars) for _ in range(1024*1024)))
len_contents = len(contents)

contents_list = contents.splitlines()

def get_opener_family_assertions(suffix):
    filename = f"{test_path}/1.{suffix}"
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
            i+=1
            assert reader._tqdm._n == i
    shell_utils.rm_rf(filename)


def test_txt():
    get_opener_family_assertions("txt")


def test_bz2():
    get_opener_family_assertions("bz2")


def test_gz():
    get_opener_family_assertions("gz")


def test_xz():
    get_opener_family_assertions("xz")


def test_lzma():
    get_opener_family_assertions("lzma")
