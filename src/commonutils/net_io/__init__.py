# raise NotImplementedError
import warnings
from typing import IO

warnings.warn("NotImplemented")


class NetIO(IO):
    pass


class NetReader(IO):
    pass


class Downloader:
    thread_number: int
    dest_filename: str
    from_url: str
    follow_30x30: bool
