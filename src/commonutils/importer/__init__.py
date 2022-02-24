"""
importer -- Wrappers of Third-party Libraries
=============================================

This module defines the importing script of third-party libraries for following circumstances:

1. The third-party libraries do not do human things in some circumstances.
For example, `tqdm` will pollute stderr if it is connected to a logger.
2. The third-party libraries is only available in some platform but not others.
For example, `pysam` can be used in GNU/Linux but not Microsoft Windows,
making it hard for people that develops under Windows get type hints.
"""
