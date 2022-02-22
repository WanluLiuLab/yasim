from abc import abstractmethod
from typing import Dict, Union

from bioutils.datastructure import Feature
from commonutils import logger
from commonutils.str_utils import to_dict

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

    def __str__(self):
        return repr(self)


class Gff3Record(BaseGtfGffRecord):
    """
    A general GTF Record.
    """

    @classmethod
    def from_string(cls, in_str: str):
        """
        To generate ONE record from existing string.

        :param in_str: Input dictionary.
        """
        global lh
        lh.debug(f'Adding {in_str}')
        line_split = in_str.split('\t')

        required_fields = line_split[0:-1]
        attributes = to_dict(line_split[-1], field_sep='=', record_sep=';', quotation_mark='\"\'', resolve_str=True)

        # Score should be an integer
        if required_fields[5] == ".":
            required_fields[5] = 0
        return Gff3Record(seqname=required_fields[0],
                          source=required_fields[1],
                          feature=required_fields[2],
                          start=int(required_fields[3]),
                          end=int(required_fields[4]),
                          score=int(float(required_fields[5])),
                          strand=(required_fields[6]),
                          frame=(required_fields[7]),
                          attribute=attributes)

    def __repr__(self):
        """
        Generate GTF string.

        :return: GTF string.
        """
        attribute_str = ""
        for k, v in self.attribute.items():
            v_str = repr(v).replace("'", '"')
            attribute_str = f"{attribute_str}{k}={v_str}; "
        return ("\t".join((
            self.seqname,
            (self.source),
            self.feature,
            str(self.start),
            str(self.end),
            str(self.score),
            self.strand,
            self.frame,
            attribute_str
        )))


class GtfRecord(BaseGtfGffRecord):
    """
    A general GTF Record.
    """

    @classmethod
    def from_string(cls, in_str: str):
        """
        To generate ONE record from existing string.

        :param in_str: Input dictionary.
        """
        global lh
        in_str = in_str.rstrip('\n\r')
        lh.debug(f'Adding {in_str}')
        line_split = in_str.split('\t')

        required_fields = line_split[0:-1]
        attributes = to_dict(line_split[-1], field_sep=' ', record_sep=';', quotation_mark='\"\'', resolve_str=True)

        # Score should be an integer
        if required_fields[5] == ".":
            required_fields[5] = 0
        return GtfRecord(seqname=required_fields[0],
                         source=required_fields[1],
                         feature=required_fields[2],
                         start=int(required_fields[3]),
                         end=int(required_fields[4]),
                         score=int(float(required_fields[5])),
                         strand=(required_fields[6]),
                         frame=(required_fields[7]),
                         attribute=attributes)

    def __repr__(self):
        """
        Generate GTF string.

        :return: GTF string.
        """
        attribute_str = ""
        for k, v in self.attribute.items():
            attribute_str = f"{attribute_str}{k} " + repr(v).replace("'", '"') + "; "
        return ("\t".join((
            self.seqname,
            (self.source),
            self.feature,
            str(self.start),
            str(self.end),
            str(self.score),
            self.strand,
            self.frame,
            attribute_str
        )))
