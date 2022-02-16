import glob
import os
from typing import List


def get_depth_from_intermediate_fasta(intermediate_fasta_dir: str) -> List[str]:
    retl = []
    for filename in glob.glob(os.path.join(intermediate_fasta_dir, "*.fa")):
        retl.append(os.path.basename(os.path.splitext(filename)[0]))
    return retl
