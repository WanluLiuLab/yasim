
from bioutils.io import GtfIterator, Gff3Iterator, get_file_type_from_suffix

file_type = get_file_type_from_suffix()
if file_type == "GTF":
    pass
elif file_type == "GFF3":
    pass
else:
    pass
