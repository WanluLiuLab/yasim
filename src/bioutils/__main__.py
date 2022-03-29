from yasim import __version__

from commonutils.libfrontend import setup_frontend

if __name__ == '__main__':
    setup_frontend(
        "bioutils.main",
        "DESC TODO",
        __version__
    )
