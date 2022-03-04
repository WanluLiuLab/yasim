from commonutils.io.tqdm_reader import get_tqdm_line_reader


def test_get_tqdm_line_reader():
    for line in get_tqdm_line_reader(__file__):
        print(line)
