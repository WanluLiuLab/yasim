"""
>>> RMSK_GFF_ATTR_REGEX.match("ID=4663330;Target=LSU-rRNA_Hsa 3753 3907").groups()
('4663330', None, 'LSU-rRNA_Hsa', '3753', '3907')
>>> RMSK_GFF_ATTR_REGEX.match("ID=129048;Target=A-rich -2 3907").groups()
('4663330', None, 'LSU-rRNA_Hsa', '3753', '3907')
>>> RMSK_GFF_ATTR_REGEX.match("ID=443;Target \\"Motif:MIR1_Amn\\" 27 207").groups()
('443', 'Motif:', 'MIR1_Amn', '27', '207')
"""

__all__ = "parse_record"

import re

from labw_utils.bioutils.parser import BaseFileIterator
from labw_utils.bioutils.record.feature import FeatureInterface, Feature
from labw_utils.commonutils.lwio.tqdm_reader import get_tqdm_line_reader
from labw_utils.typing_importer import List, Optional, Iterable, Iterator


class GFFParsingError(ValueError):
    """General GFF parsing errors."""

    pass


RMSK_GFF_ATTR_REGEX = re.compile(r"^ID=(\d+);Target[= ]\"?([^:]+:)?([^:\"]+)\"? (-?\d+) (-?\d+)$")


def parse_record(
    in_str: str,
    skip_fields: Optional[List[str]] = None,
    included_attributes: Optional[List[str]] = None,
) -> FeatureInterface:
    """
    Parse record string to :py:class:`Feature`.

    :param in_str: String to be parsed.
    :param skip_fields: Explicitly skip optional features to reduce space.
    :param included_attributes: Explicitly include attributes to reduce space. Other attributes are discarded.

    :raises GFFParsingError: On invalid record.
    """
    if skip_fields is None:
        skip_fields = []
    line_split = in_str.rstrip("\n\r").split("\t")
    if len(line_split) != 9:
        raise GFFParsingError(f"Illegal GFF record '{in_str}': Should have 9 fields, here only {len(line_split)}")
    required_fields = line_split[0:-1]
    attribute_line = line_split[-1]
    """Sample: ID=443;Target "Motif:MIR1_Amn" 27 207"""
    lm = RMSK_GFF_ATTR_REGEX.match(attribute_line)
    if lm is None:
        raise GFFParsingError(f"Cannot parse Attribute '{attribute_line}'")
    attributes = {
        "repeat_id": lm.group(1),
        "repeat_name": lm.group(3),
        "repeat_match_start": int(lm.group(4)),
        "repeat_match_end": int(lm.group(5)),
    }

    if included_attributes is not None:
        attributes = {k: v for k, v in attributes.items() if k in included_attributes}

    try:
        return Feature(
            seqname=required_fields[0],
            source=required_fields[1] if required_fields[1] != "." and "source" not in skip_fields else None,
            feature=required_fields[2] if required_fields[2] != "." and "feature" not in skip_fields else None,
            start=int(required_fields[3]),
            end=int(required_fields[4]),
            score=float(required_fields[5]) if required_fields[5] != "." and "score" not in skip_fields else None,
            strand=required_fields[6] == "+" if required_fields[6] != "." and "strand" not in skip_fields else None,
            frame=int(required_fields[7]) if required_fields[7] != "." and "frame" not in skip_fields else None,
            attribute=attributes,
        )
    except ValueError as e:
        raise GFFParsingError(f"Illegal GFF record '{in_str}'") from e


class RMSKGffIterator(BaseFileIterator, Iterable[FeatureInterface]):
    filetype: str = "RMSKGff"
    record_type = FeatureInterface

    def __iter__(self) -> Iterator[FeatureInterface]:
        for line in get_tqdm_line_reader(self.filename):
            if line.startswith("#") or line == "":
                continue
            yield parse_record(line)
