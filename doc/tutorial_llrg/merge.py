import copy
import multiprocessing
import os.path
import shutil
import subprocess
import sys
import tempfile
import uuid
from multiprocessing.managers import SyncManager
from typing import List

from labw_utils.commonutils.stdlib_helper.parallel_helper import ParallelJobExecutor


def merge(
        _d1_path: str,
        _d2_path: str,
        _o_path: str,
        _fns: List[str]
):
    if subprocess.Popen(
            [
                "pbmerge",
                "-o", _o_path,
                _d1_path,
                _d2_path
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
    ).wait() == 0:
        _fns.append(_o_path)
    else:
        print(f"MERGE {_d1_path} and {_d2_path} have errors; discarded")


if __name__ == "__main__":
    fns = sys.argv[2:]
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
                job_pool.append(multiprocessing.Process(target=merge, args=(d1_path, d2_path, o_path, fns_shared)))
            if len(this_fns) == 1:
                fns_shared.append(this_fns.pop())
            job_pool.start()
            job_pool.join()
            fns = list(fns_shared)
            print(f"Remaining: {len(fns)}")
            sm.shutdown()
            sm.join()
        shutil.move(fns[0], sys.argv[1])
