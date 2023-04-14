import glob
from yasim.helper import depth_io

if __name__ == "__main__":
    for fn in glob.glob("ce11_as_*_isoform_depth_20.tsv"):
        print(fn)
        out_fn = fn.replace("20", "25")
        depth_io.write_depth({k: v+5 for k, v in depth_io.read_depth(fn).items()}, out_fn, "TRANSCRIPT_ID")

