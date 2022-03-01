from abc import abstractmethod
from typing import Iterator, Type, Iterable, Optional, IO


class BaseIterator(Iterable):
    filename: str = ""
    filetype: str = None
    record_type: Type = None
    fd: Optional[IO]

    def __init__(self, filename: str, show_tqdm: bool = True):
        self.filename = filename
        self.show_tqdm = show_tqdm

    @abstractmethod
    def __iter__(self) -> Iterator[record_type]:
        pass

    def __repr__(self):
        return f"{self.filetype} Iterator for {self.filename}"

    __str__ = __repr__

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        return
