import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import lognorm

data = pd.read_csv("fitted_n_isoforms.csv")

count = data["count"].to_numpy()

count_extracted = []
for n in range(1, len(count) + 1):
    count_extracted.extend([n] * int(count[n - 1]))

f, loc, scale = lognorm.fit(data=count_extracted, method="MM")

count_new = []

while len(count_new) < len(count_extracted):
    n = int(lognorm.rvs(f, loc=loc, scale=scale, size=1))
    if 1 <= n < 26:
        count_new.append(n)

fig, ax = plt.subplots(1, 1)
ax.hist(count_new, bins=25, density=True, histtype="stepfilled", alpha=0.2)
ax.hist(count_extracted, bins=25, density=True, histtype="stepfilled", alpha=0.2)
plt.yscale("log")
plt.show()
print((f, loc, scale))
