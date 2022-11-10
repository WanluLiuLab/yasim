from labw_utils.commonutils.libfrontend import setup_frontend

from yasim import __version__

if __name__ == '__main__':
    setup_frontend(
        "yasim_scripts.main",
        "Supplementary scripts to yasim",
        __version__
    )
