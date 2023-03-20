import math

import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv("data.csv")
    col_names = list(df.columns)
    col_names.remove("Software")
    col_names.remove("Dataset")
    df["_sum"] = df.sum(axis=1, numeric_only=True)
    unique_software = df.Software.unique()
    unique_data = df.Dataset.unique()

    coefficient = 0.6 / len(unique_software)
    with open("coefficient.tex", "w") as writer:
        writer.write(r"\newcommand{\coefficient}{%.2f}" % coefficient)
        writer.write("\n")
    with open("table.tex", "w") as writer:
        c = "c" * (len(unique_software) + 1)
        writer.write(r"\begin{tabular}{%s}" % c)
        writer.write("\n")
        software_str = " & ".join(map(lambda _s: r"\wrapcolname{" + str.replace(_s, "_", "\_") + "}", unique_software))
        writer.write(r"    Data\textbackslash Software & %s \\" % software_str)
        writer.write("\n")
        for data in unique_data:
            this_data_str = r"    \wraprowname{%s}" % data.replace("_", "\_")
            all_software_data = df.query(f"Dataset == '{data}'")["_sum"]
            all_software_sum = sum(all_software_data)
            all_software_max = max(all_software_data)
            for software in unique_software:
                print(data, software)
                this_data_str += " & "
                this_sum = df.query(f"Dataset == '{data}' & Software == '{software}'")["_sum"].item()
                width = math.sqrt(this_sum / all_software_max)
                this_data_str += r"\wrapfig{\includegraphics[width=%.4f\linewidth]{%s}}" % (
                    width,
                    "/".join((
                        "fig",
                        f"{software}-{data}"
                    )))
            this_data_str += r"\\"
            writer.write(this_data_str)
            writer.write("\n")
        writer.write(r"\end{tabular}")
        writer.write("\n")
