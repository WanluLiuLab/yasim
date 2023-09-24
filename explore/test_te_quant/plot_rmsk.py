import os

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

CURDIR = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    df = pd.read_csv("aln/ce11_denovo_test_rmsk_prep.comp.tsv", sep="\t", quotechar="'")
    df_type = df.groupby("EVENT_TYPE").count().reset_index()
    print(df_type)
    sns.barplot(df_type, x="EVENT_TYPE", y="EVENT_ID")
    plt.savefig(os.path.join(CURDIR, "plot", "rmsk_type.png"))
    plt.cla()
