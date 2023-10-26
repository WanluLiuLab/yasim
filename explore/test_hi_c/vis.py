import cooler
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from matplotlib.ticker import EngFormatter

file = "aln/SRR9286043.100000.cool"
clr = cooler.Cooler(file)
# Define chromosome starts
chromstarts = []
for i in clr.chromnames:
    chromstarts.append(clr.extent(i)[0])

bp_formatter = EngFormatter("b")


def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_minor_formatter(bp_formatter)
    if x:
        ax.xaxis.set_minor_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis="x", rotation=45)


start, end = min(chromstarts), max(chromstarts)
mat = clr.matrix(balance=False)[start:end, start:end]

norm = LogNorm(vmin=1, vmax=mat.max())

f, ax = plt.subplots(figsize=(13, 10), nrows=1, ncols=1, sharex=False, sharey=False)

ax.set_title("Chromosomes I-V")
im = ax.matshow(
    mat,
    norm=norm,
    cmap="Reds",
)
plt.colorbar(im, ax=ax, label="Whole-genome")
ax.set_xticks(
    ticks=np.array(chromstarts) - start,
    labels=clr.chromnames,
    minor=False,
    rotation=90,
)
ax.set_yticks(
    ticks=np.array(chromstarts) - start,
    labels=clr.chromnames,
    minor=False,
    rotation=90,
)
print(clr.chromnames)
format_ticks(ax, rotate=False)

plt.tight_layout()
plt.show()
