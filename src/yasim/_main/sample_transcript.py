"""
sample_transcript.py -- Proxy to ``sample_transcript`` in ``labw_utils.bioutils``.

.. versionadded:: 3.1.4

.. versiondeprecated:: 3.1.5
    Use that in ``labw_utils.bioutils`` instead.
"""
from labw_utils.commonutils.stdlib_helper.logger_helper import get_logger
from labw_utils.typing_importer import List

try:
    from labw_utils.bioutils.main.sample_transcript import main as bc_main
except ImportError:
    from labw_utils.bioutils._main.sample_transcript import main as bc_main

try:
    from labw_utils.bioutils.main.transcribe import create_parser as create_parser
except ImportError:
    try:
        from labw_utils.bioutils._main.transcribe import create_parser as create_parser
    except ImportError:
        create_parser = None

_lh = get_logger(__name__)


def main(args: List[str]):
    _lh.warning(
        "DEPRECATION WARNING: "
        "`sample_transcript` at `yasim` had been deprecated. Use those in `labw_utils.bioutils` instead"
    )
    bc_main(args)
