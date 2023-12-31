import sys
from typing import List


def main(_: List[str]):
    name_cache = None
    for line in sys.stdin:
        if name_cache is None:
            name_cache = line
        else:
            if len(line) >= 256:
                try:
                    sys.stdout.write(name_cache)
                    sys.stdout.write(line)
                except BrokenPipeError:
                    break
            name_cache = None
