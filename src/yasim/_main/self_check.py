"""
self_check.py -- Check whether YASIM installation is complete.

.. versionadded:: 3.1.6
"""
import importlib
import shutil

from labw_utils.typing_importer import List, Callable, Any


def get_version(
    pkg_name: str,
    version_hook: Callable[[Any], str] = lambda m: getattr(m, "__version__"),
):
    try:
        mod = importlib.import_module(pkg_name)
        ver = version_hook(mod)
    except (
        ImportError,
        ValueError,
        AttributeError,
        RuntimeError,
        SystemError,
        TypeError,
    ):
        ver = "ERR"
    print(f"{pkg_name}: {ver}")


def get_exec(exec_name: str):
    exec_path = shutil.which(exec_name)
    if exec_path is None:
        exec_path = "ERR"
    print(f"{exec_name}: {exec_path}")


def main(args: List[str]):
    _ = args
    del args

    print(f"YASIM version info:")
    for pkg_name in [
        "labw_utils",
        "yasim",
        "pandas",
        "numpy",
        "scipy",
        "tqdm",
        "joblib",
        "pysam",
        "jinja2",
        "h5py",
        "setuptools",
        "seaborn",
    ]:
        get_version(pkg_name)
    get_version("matplotlib", lambda m: m._version.version)

    for exec_name in [
        "badread",
        "pbsim",
        "pbsim2",
        "pbsim3",
        "dwgsim",
        "art_illumina",
        "ccs",
        "samtools",
    ]:
        get_exec(exec_name)
