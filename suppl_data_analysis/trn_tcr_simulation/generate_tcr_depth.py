from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader
from yasim.helper.depth import write_depth


if __name__ == "__main__":
    barcode_path = "barcode.txt"
    depth_db = {}
    for barcode in get_tqdm_line_reader(barcode_path):
        depth_db[f"{barcode}:A"] = 500
        depth_db[f"{barcode}:B"] = 500

    write_depth(depth_db, "tcr_depth.tsv")
