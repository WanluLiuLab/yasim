# ==============================================================================
#  Copyright (C) 2021. tetgs authors
#
#  This file is a part of tetgs, which is licensed under MIT,
#  a copy of which can be obtained at <https://opensource.org/licenses/MIT>.
#
#  NAME: setup.py -- Installer
#
#  VERSION HISTORY:
#  2021-09-11 0.1  : Purposed and added by YU Zhejian.
#
# ==============================================================================
import glob
import os.path
import sys

import setuptools
from setuptools import setup

PKG_NAME = "yasim"

ROOT_DIR = os.path.dirname(__file__)
sys.path.append(os.path.join(ROOT_DIR, "src"))

install_requires = []

with  open('requirements.txt', 'rt', encoding='utf-8') as reader:
    for line in reader:
        if not line.startswith('#'):
            install_requires.append(line.strip())

from yasim import __version__

with  open('Readme.md', 'rt', encoding='utf-8') as reader:
    long_description = reader.read()

setup(
    name=PKG_NAME,
    version=__version__,
    author="YU Zhejian",
    author_email="Zhejian.19@intl.zju.edu.cn",
    description=f"{PKG_NAME} -- A Simulator for Alternative Splicing and Differentially Expressed Gene",
    long_description=long_description,
    long_description_content_type='text/markdown',
    url="https://github.com/pypa/sampleproject",  # TODO
    project_urls={  # TODO
        "Bug Tracker": "",  # TODO
        "Documentations": ""  # TODO
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX",
        "Programming Language :: R"
    ],
    python_requires=">=3.7",
    packages=setuptools.find_packages(
        where='src',
        include=['*'],
    ),
    package_dir={"": 'src'},
    package_data={
        '': glob.glob(os.path.join(ROOT_DIR, "src", "yasim", "llrg_adapter", "**"), recursive=True),
    },
    install_requires=install_requires
    # FIXME: Errors when adding to sdists.
)
