raise NotImplementedError


class NetIO(IO):
    pass


class NetReader(IO):
    pass


class Downloader:
    thread_number: int
    dest_filename: str
    from_url: str
    follow_30x30: bool
