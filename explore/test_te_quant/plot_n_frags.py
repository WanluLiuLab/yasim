import os

import seaborn as sns
import matplotlib.pyplot as plt

from yasim.helper.translation_instruction import TranslationInstruction

CURDIR = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    ti = TranslationInstruction.from_json(os.path.join(CURDIR, "sim", "ce11_denovo_test.json"))
    retl = []
    for t in ti.transcripts.values():
        retl.append(len(t.l))
    sns.histplot(x=retl)
    plt.savefig(os.path.join(CURDIR, "plot", "n_frags.png"))
