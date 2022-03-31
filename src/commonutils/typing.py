# import mmap
import io
from io import TextIOWrapper, StringIO, BufferedIOBase
from os import PathLike
from typing import Union, IO

PathType = Union[PathLike, str]

FDType = Union[IO, io.IOBase, StringIO]

PathOrFDType = Union[PathType, FDType]



