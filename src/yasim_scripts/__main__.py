"""
TODO docs

.. versionadded:: 3.1.5
"""

from labw_utils.commonutils.libfrontend import setup_frontend

from yasim import __version__
from yasim_scripts import description

if __name__ == "__main__":
    setup_frontend(f"{__package__}._main", description, __version__)
