"""
The Low-Level Parser of GTF and GFF3 Files
==========================================

This module includes record iterators of GTF and GFF3 files,
and a tree parser of GFF3 files.

The parsers are organized in following ways:

The "record iterators" are like SAX parsers when paring XML files.
It parses GTF/GFF3 files into a one-direction iterator of :py:class:`GtfRecord` or :py:class:`Gff3Record`.
The memory cost of record iterator is small.

Since GFF3 files are organized in a Child-Parent way,
a DOM-like GFF3 tree (actually Graph) parser is implemented.
This parser parses GFF3 files into a Directed-Acyclic Graph (DAG).

Warnings:
---------

This module does not:

* Parse a GTF/GFF3 record string into :py:class:`GtfRecord`.
See :py:class:`bioutils.typing.feature.GtfRecord` for this feature.

* Parse a GTF/GFF3 file into a three-tier Exon-Transcript-GeneView structure.
See :py:class:`bioutils.datastructure.gene_view.GeneView` for this feature.
"""

from collections import defaultdict
from typing import Dict, Iterable, Iterator, Union, Optional, List, TextIO

from bioutils.typing.feature import Feature, GFF3_TOPLEVEL_NAME, Gff3Record, GtfRecord
from bioutils.typing.misc import BaseIterator
from commonutils.importer.tqdm_importer import tqdm
from commonutils.io.safe_io import get_writer
from commonutils.io.tqdm_reader import get_tqdm_line_reader


class GtfIterator(BaseIterator):
    filetype: str = "GTF"
    record_type = GtfRecord

    def __iter__(self) -> Iterator[GtfRecord]:
        for line in get_tqdm_line_reader(self.filename):
            if line.startswith('#') or line == '':
                continue
            yield GtfRecord.from_string(line)


class Gff3Iterator(BaseIterator):
    filetype: str = "GFF3"
    record_type = Gff3Record

    def __iter__(self) -> Iterator[Gff3Record]:
        for line in get_tqdm_line_reader(self.filename):
            if line.startswith('##FASTA'):
                return
            if line.startswith('#') or line == '':
                continue
            yield Gff3Record.from_string(line)


class _FeatureWriter:
    df: TextIO

    @staticmethod
    def write_iterator(
            iterable: Union[GtfIterator, GtfIterator, Iterator[Feature]],
            output_filename: str,
            prefix_annotations: Optional[List[str]] = None
    ):
        with _FeatureWriter(output_filename) as writer:
            if prefix_annotations is not None:
                for annotation in prefix_annotations:
                    writer.write_comment(annotation)
            for feature in iterable:
                writer.write_feature(feature)

    def __init__(self, output_filename: str):
        self.output_filename = output_filename
        self.fd = get_writer(self.output_filename)

    def write_feature(self, feature: Feature, **kwargs):
        self.fd.write(feature.format_string(**kwargs) + "\n")

    def write_comment(self, comment: str):
        self.fd.write('#' + comment + "\n")

    def close(self):
        self.fd.close()

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        return

    def tell(self) -> int:
        try:
            return self.fd.tell()
        except OSError:
            return -1

    def __repr__(self) -> str:
        return f"Feature writer of {self.output_filename} at {self.tell()}"

    __str__ = __repr__


class GtfWriter(_FeatureWriter):
    pass


class Gff3Writer(_FeatureWriter):

    def __init__(self, output_filename: str):
        super().__init__(output_filename=output_filename)
        self.write_comment("#gff-version 3")  # This comment should be wrieen at the first line.


class Gff3Tree:
    filename: str

    _flat_record: Dict[str, Gff3Record]
    """An ID -> Record mapping"""

    _id_tree: Dict[str, List[str]]
    """An parent ID -> child ID mapping"""

    def __init__(self, filename: str):
        self.filename = filename
        self._flat_record = {}
        self._id_tree = defaultdict(lambda: [])
        self._load()
        self._parse_tree()

    def _load(self):
        self._flat_record = {}
        for gff3_record in Gff3Iterator(self.filename):
            self._flat_record[gff3_record.id] = gff3_record

    def _parse_tree(self):
        for k, v in tqdm(iterable=self._flat_record.items(), desc="Parsing GFF3 to a Tree..."):
            self._id_tree[v.parent_id].append(k)

    def get(self, gff3_id: str) -> Gff3Record:
        if gff3_id == GFF3_TOPLEVEL_NAME:
            raise ValueError(f"Cannot reach for {GFF3_TOPLEVEL_NAME} -- It is toplevel virtual ID")
        return self._flat_record[gff3_id]

    def get_child_ids(self, gff3_id: str) -> Iterable[str]:
        return self._id_tree[gff3_id]

    def get_all_ids(self) -> Iterable[str]:
        return self._flat_record.keys()

    def get_toplevel_ids(self) -> Iterable[str]:
        return self.get_child_ids(GFF3_TOPLEVEL_NAME)

    def __len__(self) -> int:
        return len(self._flat_record)

    def get_parent_id(self, gff3_id: str) -> str:
        return self.get(gff3_id).parent_id

    def get_backtrack_feature(self, gff3_id: str, parent_feature_name: str) -> str:
        while True:
            gff3_id = self.get_parent_id(gff3_id)
            if self.get(gff3_id).feature == parent_feature_name:
                return gff3_id

        raise ValueError("Backtrack failed!")
