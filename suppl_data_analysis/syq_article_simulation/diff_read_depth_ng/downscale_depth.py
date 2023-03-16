import sys

from yasim.helper.depth import DepthType, read_depth, write_depth

if __name__ == "__main__":
    d = read_depth(sys.argv[1])
    d_mutated = {
        k: d[k] * float(sys.argv[3])
        for k in d.keys()
    }
    d_mutated_filtered = {
        k: 0 if v < 0.01 else v
        for k, v in d_mutated.items()
    }
    write_depth(d_mutated_filtered, sys.argv[2])
