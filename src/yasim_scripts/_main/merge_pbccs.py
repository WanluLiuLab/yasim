"""
merge_pbccs.py -- Merge BAMs created by pbccs.
"""

import argparse
import copy
import glob
import multiprocessing
import os.path
import shutil
import tempfile
import uuid
from multiprocessing.managers import SyncManager
from typing import List

from labw_utils.commonutils.stdlib_helper.argparse_helper import (
    ArgumentParserWithEnhancedFormatHelp,
)
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.commonutils.stdlib_helper.parallel_helper import (
    ParallelJobExecutor,
    easyexec,
)

__all__ = (
    "create_parser",
    "main",
)

_lh = get_logger(__name__)


def create_parser() -> argparse.ArgumentParser:
    parser = ArgumentParserWithEnhancedFormatHelp(
        prog="python -m yasim_scripts merge_pbccs",
        description=__doc__.splitlines()[1],
    )
    parser.add_argument(
        "-o",
        "--out",
        help="Output BAM file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-e",
        "--pbmerge_path",
        help="Path to pbmerge",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--input_bam_glob",
        type=str,
        required=True,
        help="Glob expression for BAM files that will be merged.",
    )
    return parser


def merge(_d1_path: str, _d2_path: str, _o_path: str, _fns: List[str]) -> None:
    if (
        easyexec(["pbmerge", "-o", _o_path, _d1_path, _d2_path], raise_on_error=False)
        == 0
    ):
        _fns.append(_o_path)
    else:
        _lh.error(f"MERGE %s and %s have errors; discarded", _d1_path, _d2_path)


def main(args: List[str]) -> int:
    argv = create_parser().parse_args(args)
    pbmerge_path = argv.pbmerge_path
    if pbmerge_path is None:
        pbmerge_path = shutil.which("pbmerge")
    if pbmerge_path is None:
        _lh.error("pbmerge not found!")
        return 1
    fns = list(glob.glob(argv.input_bam_glob))
    i = 0
    with tempfile.TemporaryDirectory() as tmpdir:
        while len(fns) > 1:
            sm = SyncManager()
            sm.start()
            i += 1
            print(f"Iter {i} with {len(fns)} files remaining")
            this_fns: List[str] = copy.deepcopy(fns)
            fns_shared = sm.list()
            job_pool = ParallelJobExecutor()
            while len(this_fns) >= 2:
                d1_path, d2_path = this_fns.pop(), this_fns.pop()
                o_path = os.path.join(tmpdir, str(uuid.uuid4()) + ".bam")
                job_pool.append(
                    multiprocessing.Process(
                        target=merge, args=(d1_path, d2_path, o_path, fns_shared)
                    )
                )
            if len(this_fns) == 1:
                fns_shared.append(this_fns.pop())
            job_pool.start()
            job_pool.join()
            fns = list(fns_shared)
            print(f"Remaining: {len(fns)}")
            sm.shutdown()
            sm.join()
        shutil.move(fns[0], argv.out)
    return 0
