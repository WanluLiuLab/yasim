import os

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

CURDIR = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    df = pd.read_csv(os.path.join(CURDIR, "assmb_final", "transrate", "assemblies.csv"))
    df["assembly"] = list(map(lambda p: os.path.basename(p), df["assembly"]))
    plt.figure(figsize=(10, 10))

    assemb_ns_y_values = list(map(lambda x: "n" + str(x), range(10, 100, 20)))
    assemb_ns = df.loc[:, ["assembly", *assemb_ns_y_values]].melt(
        id_vars="assembly", value_vars=assemb_ns_y_values, value_name="val", var_name="cat"
    )
    sns.lineplot(data=assemb_ns, x="cat", y="val", hue="assembly")
    plt.savefig(os.path.join(CURDIR, "plot", "assemb_ns.png"))
    plt.cla()

    assemb_cov_values = ["p_cov25", "p_cov50", "p_cov75", "p_cov85", "p_cov95", "reference_coverage"]
    assemb_cov = df.loc[:, ["assembly", *assemb_cov_values]].melt(
        id_vars="assembly", value_vars=assemb_cov_values, value_name="val", var_name="cat"
    )
    sns.lineplot(data=assemb_cov, x="cat", y="val", hue="assembly")
    plt.savefig(os.path.join(CURDIR, "plot", "assemb_cov.png"))
    plt.cla()

    assemb_crbb_values = ["p_contigs_with_CRBB", "p_refs_with_CRBB"]
    assemb_crbb = df.loc[:, ["assembly", *assemb_crbb_values]].melt(
        id_vars="assembly", value_vars=assemb_crbb_values, value_name="val", var_name="cat"
    )
    assemb_crbb_plot = sns.barplot(data=assemb_crbb, x="assembly", y="val", hue="cat")
    assemb_crbb_plot.set_xticklabels(
        assemb_crbb_plot.get_xticklabels(),
        rotation=30,
        horizontalalignment="right",
    )
    plt.savefig(os.path.join(CURDIR, "plot", "assemb_crbb.png"))
    plt.cla()

    assemb_crbb_plot = sns.barplot(data=df, x="assembly", y="CRBB_hits")
    assemb_crbb_plot.set_xticklabels(
        assemb_crbb_plot.get_xticklabels(),
        rotation=30,
        horizontalalignment="right",
    )
    plt.savefig(os.path.join(CURDIR, "plot", "assemb_crbb_n.png"))
    plt.cla()
