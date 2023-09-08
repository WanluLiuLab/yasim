"""
maf_parser.py -- Parser for LAST aligned/PBSIM MAF.

.. versionadded:: 3.1.5
"""
__all__ = ("MAF_RECORD_REGEX", "MafRecordType", "maf_parse")

import re

from labw_utils.commonutils.lwio.safe_io import get_reader
from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import Tuple, Iterable

MafRecordType = Tuple[str, str, str, str]
"""
Name1, Name2, Alignment1, Alignment2

.. versionadded:: 3.1.5
"""

MAF_RECORD_REGEX = re.compile(r"^s +(\S+) +(\d+) +(\d+) +([+-]) +(\d+) +(\S+)$")
"""
Regex to MAF record

.. versionadded:: 3.1.5
"""

_lh = get_logger(__name__)


def maf_parse(maf_path: str, show_tqdm: bool = True) -> Iterable[MafRecordType]:
    """
    Parse MAF into iterable of ``MafRecordType``.

    :param maf_path: Path to input MAF file.
    :param show_tqdm: Whether to show progress bar.

    .. versionadded:: 3.1.5
    """
    num_error = 0
    num_record = 0
    if show_tqdm:
        rf = get_tqdm_line_reader
    else:
        rf = get_reader
    with rf(maf_path) as reader:
        while True:
            line1 = reader.readline()
            if line1 == "":
                break
            if not line1.startswith("s"):
                continue
            else:
                line2 = reader.readline().rstrip("\n")
                line1 = line1.rstrip("\n")
            lm1 = MAF_RECORD_REGEX.match(line1)
            lm2 = MAF_RECORD_REGEX.match(line2)
            if lm1 is None or lm2 is None:
                num_error += 1
                continue
            g1 = lm1.groups()
            g2 = lm2.groups()
            yield g1[0], g2[0], g1[5], g2[5]
            num_record += 1
    _lh.debug("Finished with %d errors and %d records", num_error, num_record)
