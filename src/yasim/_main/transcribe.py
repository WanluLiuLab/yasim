"""
transcribe.py -- Proxy to ``transcribe`` in ``labw_utils.bioutils``.
"""
from typing import List

from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger

try:
    from labw_utils.bioutils.main.transcribe import main as bc_main
except ImportError:
    from labw_utils.bioutils._main.transcribe import main as bc_main

try:
    from labw_utils.bioutils.main.transcribe import create_parser as create_parser
except ImportError:
    try:
        from labw_utils.bioutils._main.transcribe import create_parser as create_parser
    except ImportError:
        create_parser = None

_lh = get_logger(__name__)


def main(args: List[str]):
    _lh.warning("DEPRECATION WARNING: "
                "`transcribe` at `yasim` had been deprecated. Use those in `labw_utils.bioutils` instead")
    bc_main(args)
