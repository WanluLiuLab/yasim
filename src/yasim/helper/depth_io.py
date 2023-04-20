"""
depth_io.py -- GEP Datastructure and Utils

Here contains Datastructure and I/O utilities for Gene Expression Profile (GEP).
"""

__all__ = (
    "DepthType",
    "write_depth",
    "read_depth",
    "DepthParsingException"
)

from labw_utils.typing_importer import Dict, Literal

from labw_utils.commonutils.importer.tqdm_importer import tqdm
from labw_utils.commonutils.io.safe_io import get_writer, get_reader
from labw_utils.commonutils.io.tqdm_reader import get_tqdm_line_reader

DepthType = Dict[str, float]
"""DGE type, is transcript_id -> coverage"""


class DepthParsingException(RuntimeError):
    ...


def write_depth(
        depth_data: DepthType,
        dst_depth_file_path: str,
        feature_name: Literal["GENE_ID", "TRANSCRIPT_ID"],
        show_tqdm: bool = True
):
    """
    Write Depth information to file

    :param depth_data: Data in depth.
    :param dst_depth_file_path: Output TSV path. Can be compressed.
    :param feature_name: Name of deature. Should be ``GENE_ID`` or ``TRANSCRIPT_ID``.
    :param show_tqdm: Whether to show a progress bar.
    """
    with get_writer(dst_depth_file_path) as writer:
        writer.write(f"{feature_name}\tDEPTH\n")
        it = depth_data.items()
        if show_tqdm:
            it = tqdm(it, desc=f"Writing to {dst_depth_file_path}")
        for transcript_id, d in it:
            writer.write(f"{transcript_id}\t{d}\n")


def read_depth(
        src_depth_file_path: str,
        show_tqdm: bool = True
) -> DepthType:
    """
    Read depth information from a file.

    :param src_depth_file_path: Path to source depth file. Can be compressed.
    :param show_tqdm: Whether to show a progress bar.
    :return: Retrived depth information
    """
    try:
        retd = {}
        if show_tqdm:
            reader = get_tqdm_line_reader(src_depth_file_path)
        else:
            reader = get_reader(src_depth_file_path)
        _ = next(reader)  # Skip line 1
        for line in reader:
            line = line.strip()
            lkv = line.split("\t")
            retd[lkv[0]] = float(lkv[1])
        reader.close()
        return retd
    except (OSError, KeyError, IndexError) as e:
        raise DepthParsingException(f"Cannot parse {src_depth_file_path}") from e
