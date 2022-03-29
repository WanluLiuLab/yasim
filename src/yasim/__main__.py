
from yasim import \
    __version__ as yasim_ver, \
    description

from commonutils.libfrontend import setup_frontend

if __name__ == '__main__':
    setup_frontend(
        "yasim.main",
        description,
        yasim_ver,
    )
