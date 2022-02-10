from typing import Dict, Union, Any

from bioutils.datastructure import Feature
from bioutils.datastructure._base_gff_gtf_record import BaseGtfGffRecord
from commonutils import logger
from commonutils.str_utils import to_dict
from bioutils.io.gtf._gtf_attribute_parser import parse as parse_gtf_attrs

lh = logger.get_logger(__name__)

__version__ = 0.1

GTFAttributeType = Dict[str, Union[str, int, float, bool, None]]

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
