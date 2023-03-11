from typing import List

try:
    from labw_utils.bioutils.main.transcribe import main as bc_main
except ImportError:
    from labw_utils.bioutils._main.transcribe import main as bc_main


def main(args: List[str]):
    bc_main(args)
