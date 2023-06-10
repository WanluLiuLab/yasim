import sys

import numpy as np

from labw_utils.commonutils.libfrontend import setup_basic_logger
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.mlutils.ndarray_helper import describe
from yasim.helper.depth_io import write_depth, read_depth

_lh = get_logger(__name__)

if __name__ == "__main__":
    setup_basic_logger()
    d = read_depth(sys.argv[1])
    d_mutated = {
        k: max(1, d[k] * float(sys.argv[3]))
        for k in d.keys()
    }
    _lh.info(
        "Generation of isoform-level depth: Final distribution: %s",
        describe(np.array(list(d_mutated.values())))
    )
    write_depth(d_mutated, sys.argv[2], "TRANSCRIPT_ID")
