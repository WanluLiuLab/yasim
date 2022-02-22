from commonutils.tqdm_utils import tqdm_line_iterator


def test_tqdm_line_iterator():
    for line in tqdm_line_iterator(__file__):
        print(line)
