from abc import abstractmethod
from typing import Dict, Union, Any

from bioutils.datastructure import Feature
from commonutils import logger

lh = logger.get_logger(__name__)

__version__ = 0.1

GTFAttributeType = Dict[str, Union[str, int, float, bool, None]]


class BaseGtfGffRecord(Feature):
    """
    A general GTF Record.
    """

    __slots__ = ('attribute',)
    attribute: GTFAttributeType

    def __init__(self,
                 seqname: str,
                 source: str,
                 feature: str,
                 start: int,
                 end: int,
                 score: float,
                 strand: str,
                 frame: str,
                 attribute: GTFAttributeType):
        """
        The filenames are named after Ensembl specifications.

        .. warning::
            Ensembl uses different way to represent 5'UTR.
        """
        super(BaseGtfGffRecord, self).__init__(
            seqname=seqname,
            source=source,
            feature=feature,
            start=start,
            end=end,
            score=score,
            strand=strand,
            frame=frame
        )
        self.attribute = attribute

    @classmethod
    def from_string(cls, in_str: str):
        """
        To generate ONE record from existing string.

        :param in_str: Input dictionary.
        """
        pass

    @abstractmethod
    def __repr__(self):
        pass

    __str__ = __repr__
