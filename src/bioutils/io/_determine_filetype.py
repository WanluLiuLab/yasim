# Common file suffixes.
common_suffixes = {
    "GTF": ('.gtf', '.gff'),
    "GFF3": ('.gff3',),
    "BED": ('.bed',),
    "RMSK": ('rmsk.txt',),
    "FASTA": ('.fasta', '.fa'),
}

archive_suffixes = {
    "GZ": (".gz", ".gzip"),
    "LZMA": (".xz", ".lzma"),
    "LZ4": (".lz4",),
    "LZIP": (".lz",),
    "LZOP": (".lzop"),
    "BZ2": (".bz2",),
    "COMPRESS": (".z",),
    "BROTLI": (".brotli", ".br"),
    "ZSTD": (".zst", ".zstd"),
    "ZIP": (".zip",),
    "RAR": (".rar", ".rar5"),
    "7Z": (".7z",)
}


def get_file_type_from_suffix(filename: str) -> str:
    filename = filename.lower()
    archive_suffix_all_removed = False
    while not archive_suffix_all_removed:
        archive_suffix_all_removed = True
        for real_suffixes in archive_suffixes.values():
            for real_suffix in real_suffixes:
                if filename.endswith(real_suffix):
                    filename = filename.rstrip(real_suffix)
                    archive_suffix_all_removed = False

    for standard_suffix, real_suffixes in common_suffixes.items():
        for real_suffix in real_suffixes:
            if filename.endswith(real_suffix):
                return standard_suffix
    return 'UNKNOWN'
