# import mmap
from io import TextIOWrapper, StringIO, BufferedIOBase
from os import PathLike
from typing import Union, IO, AnyStr

PathType = Union[PathLike, str]

FDType = Union[BufferedIOBase, TextIOWrapper, IO]

PathOrFDType = Union[PathType, FDType, StringIO]
