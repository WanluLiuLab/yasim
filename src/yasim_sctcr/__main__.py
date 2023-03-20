from labw_utils.commonutils.libfrontend import setup_frontend

from yasim import __version__

if __name__ == '__main__':
    setup_frontend(
        "yasim_sctcr._main",
        "YASIM for Single-Cell TCR-Seq Simulation",
        __version__
    )
