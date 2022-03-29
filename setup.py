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

from yasim import \
    __version__ as yasim_ver, \
        __author__ as yasim_auth, \
            author_email as yasim_auth_email
                

with  open('Readme.md', 'rt', encoding='utf-8') as reader:
    long_description = reader.read()

setup(
    name=PKG_NAME,
    version=yasim_ver,
    author=yasim_auth,
    author_email=yasim_auth_email,
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
    install_requires=install_requires,
    scripts=list(glob.glob(os.path.join(ROOT_DIR, "src", "bin", "*")))
    # FIXME: Errors when adding to sdists.
)
