import copy
import multiprocessing
import os.path
import shutil
import sys
import tempfile
import uuid
from typing import List

from pbcore.io import DataSet


def merge(d1_path: str, d2_path: str, o_path: str):
    (
        DataSet(d1_path).
        merge(
            DataSet(d2_path)
        ).
        write(o_path)
    )


if __name__ == "__main__":
    fns = sys.argv[2:]
    i = 0
    with tempfile.TemporaryDirectory() as tmpdir:
        while len(fns) > 1:
            i += 1
            print(f"Iter {i} with {len(fns)} files remaining")
            this_fns:List[str] = copy.deepcopy(fns)
            fns = []
            job_pool = []
            while len(this_fns) >= 2:
                d1_path, d2_path = this_fns.pop(), this_fns.pop()
                o_path = os.path.join(tmpdir, str(uuid.uuid4()) + ".xml")
                fns.append(o_path)
                job_pool.append(multiprocessing.Process(target=merge, args=(d1_path, d2_path, o_path)))
                job_pool[-1].start()
            if len(this_fns) == 1:
                fns.append(this_fns.pop())
            for job in job_pool:
                job.join()
        shutil.move(fns[0], sys.argv[1])
