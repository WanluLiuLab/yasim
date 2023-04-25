import argparse
import logging
import math
import os
import shutil
import subprocess
import sys

import pandas as pd

FILE_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

_s_handler = logging.StreamHandler()
_s_handler.setLevel("INFO")

logging.basicConfig(
    level="INFO",
    handlers=[
        _s_handler
    ]
)

_lh = logging.getLogger("latex_pie_chart_gen")


def gen_figs(
        src_data_csv_file_path: str,
        dst_fig_dir_path: str,
        plot_log_file_path: str,
        tmp_dir_path: str,
        r_path: str
) -> None:
    """
    Generate Figrues using R
    """
    _lh.info("GEN_FIGS: START")
    cmd = [
        r_path,
        os.path.join(FILE_DIR_PATH, "plot.R"),
        "--src_data_csv_file_path", src_data_csv_file_path,
        "--dst_fig_dir_path", dst_fig_dir_path
    ]
    _lh.info("CMD: %s", " ".join(cmd))
    with open(plot_log_file_path, "ab") as plot_log_file_writer:
        retv = subprocess.Popen(
            cmd,
            cwd=tmp_dir_path,
            stdout=plot_log_file_writer,
            stderr=plot_log_file_writer,
            stdin=subprocess.DEVNULL
        ).wait()
        if retv != 0:
            _lh.error(
                "Failed in creating figures (Return %d). See %s for details.",
                retv,
                plot_log_file_path
            )
            sys.exit(1)
    _lh.info("GEN_FIGS: FIN")


def gen_table(
        src_data_csv_file_path: str,
        src_fig_dir_path: str,
        dst_coeff_tex_file_path: str,
        dst_table_tex_file_path: str,
) -> None:
    """
    Generate Table of Figures in LaTeX.

    :return: None
    """
    _lh.info("GEN_TABLE: START")
    df = pd.read_csv(src_data_csv_file_path)
    col_names = list(df.columns)
    col_names.remove("Software")
    col_names.remove("Dataset")
    df["_sum"] = df.sum(axis=1, numeric_only=True)
    unique_software = df.Software.unique()
    unique_data = df.Dataset.unique()

    coefficient = 0.6 / len(unique_software)
    with open(dst_coeff_tex_file_path, "w") as writer:
        writer.write(r"\newcommand{\coefficient}{%.2f}" % coefficient)
        writer.write("\n")
    with open(dst_table_tex_file_path, "w") as writer:
        c = "c" * (len(unique_software) + 1)
        writer.write(r"\begin{tabular}{%s}" % c)
        writer.write("\n")
        software_str = " & ".join(
            map(lambda _s: r"\wrapcolname{" + str.replace(_s, "_", "\_") + "}", unique_software)
        )
        writer.write(r"    Data\textbackslash Software & %s \\" % software_str)
        writer.write("\n")
        for data in unique_data:
            this_data_str = r"    \wraprowname{%s}" % data.replace("_", "\_")
            all_software_data = df.query(f"Dataset == '{data}'")["_sum"]
            all_software_max = max(all_software_data)
            for software in unique_software:
                _lh.debug("GEN_TABLE: Parsing %s, %s", data, software)
                this_data_str += " & "
                try:
                    this_sum = df.query(f"Dataset == '{data}' & Software == '{software}'")["_sum"].item()
                except ValueError:
                    _lh.error("Software %s data %s have !=1 data!", software, data)
                    sys.exit(1)
                width = math.sqrt(this_sum / all_software_max)
                # LaTeX does not support pathsep in \\
                fig_file_path = r"/".join((src_fig_dir_path, f"{software}-{data}.pdf"))
                if not os.path.exists(fig_file_path):
                    raise FileNotFoundError(f"File {fig_file_path} not found!")
                this_data_str += r"\wrapfig{\includegraphics[width=%.4f\linewidth]{%s}}" % (width, fig_file_path)
            this_data_str += r"\\"
            writer.write(this_data_str)
            writer.write("\n")
        writer.write(r"\end{tabular}")
        writer.write("\n")
        _lh.info("GEN_TABLE: FIN")


def compile_latex(
        tmp_dir_path: str,
        latex_path: str,
        compile_latex_log_file_path: str
) -> None:
    _lh.info("COMPILE_LATEX: START")
    with open(compile_latex_log_file_path, "ab") as plot_log_file_writer:
        retv = subprocess.Popen(
            [
                latex_path,
                "main"
            ],
            cwd=tmp_dir_path,
            stdout=plot_log_file_writer,
            stderr=plot_log_file_writer,
            stdin=subprocess.DEVNULL
        ).wait()
        if retv != 0:
            _lh.error(
                "Failed in creating PDF (Return %d). See %s for details.",
                retv,
                compile_latex_log_file_path
            )
            sys.exit(1)
    _lh.info("COMPILE_LATEX: FIN")


if __name__ == "__main__":
    _lh.info("latex_pie_chart_gen -- Generate pie chart table")
    _lh.info("args: %s", " ".join(sys.argv))
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--data_csv_file_path")
    parser.add_argument("-o", "--out_pdf_file_path")
    parser.add_argument("-t", "--tmp_dir_path", default=os.path.abspath("tmp"))
    args = parser.parse_args(sys.argv[1:])

    _lh.info("Preparing...")
    data_csv_file_path = os.path.abspath(args.data_csv_file_path)
    _tmp_dir_path = os.path.abspath(args.tmp_dir_path)
    fig_dir_path = os.path.join(_tmp_dir_path, "figs")
    os.makedirs(_tmp_dir_path, exist_ok=True)

    _r_path = shutil.which("Rscript")
    if _r_path is None:
        _lh.error(
            "Executable of Rscript not found!"
        )
        sys.exit(1)
    _latex_path = shutil.which("pdflatex")
    if _latex_path is None:
        _lh.error(
            "Executable of pdflatex not found!"
        )
        sys.exit(1)

    shutil.copy(
        os.path.join(FILE_DIR_PATH, "main.tex"),
        os.path.join(_tmp_dir_path, "main.tex")
    )
    gen_figs(
        src_data_csv_file_path=data_csv_file_path,
        dst_fig_dir_path=fig_dir_path,
        plot_log_file_path=os.path.join(_tmp_dir_path, "gen_figs.log"),
        tmp_dir_path=_tmp_dir_path,
        r_path=_r_path

    )
    gen_table(
        src_data_csv_file_path=data_csv_file_path,
        src_fig_dir_path=fig_dir_path,
        dst_coeff_tex_file_path=os.path.join(_tmp_dir_path, "coefficient.tex"),
        dst_table_tex_file_path=os.path.join(_tmp_dir_path, "table.tex")
    )
    compile_latex(
        tmp_dir_path=_tmp_dir_path,
        latex_path=_latex_path,
        compile_latex_log_file_path=os.path.join(_tmp_dir_path, "compile_latex.log")
    )
    shutil.move(os.path.join(_tmp_dir_path, "main.pdf"), args.out_pdf_file_path)
    _lh.info("Finished")
