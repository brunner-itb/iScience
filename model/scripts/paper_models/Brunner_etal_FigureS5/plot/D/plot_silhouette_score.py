import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
import seaborn as sns
from sklearn.metrics import silhouette_samples, calinski_harabasz_score

from thesis.scripts.paper_models.utilities.plotting_rc import rc_ticks
from thesis.cellBehaviourUtilities.bridson_sampling import bridson
from thesis.cellBehaviourUtilities.grid_clustering import make_clusters

s = 200
# cell_grid_positions = get_cell_grid_positions(s, s, s)

cell_grid_positions = bridson(30,[0,0,0],[s,s,s], density_function=lambda x: 12)

apcs = np.array([
    np.mean(cell_grid_positions,axis = 0),
])

apcs = np.array([
    [130,130,130],
    [260,260,260]

])
fractions_dict = {"Th": 0.9, "Tsec": 0.1}
assert np.sum(list(fractions_dict.values())) <= 1

x = np.linspace(0,1,25)
# x = [1]
y = []

for rep in np.arange(20):
    print("replicat:", rep)
    for cs in x:
        cluster_strengths = [0,cs]

        cell_type = make_clusters(cell_grid_positions, apcs, fractions_dict, cluster_strengths)

        df = pd.DataFrame(
            {"x": cell_grid_positions[:, 0], "y": cell_grid_positions[:, 1], "z": cell_grid_positions[:, 2],
             "type": cell_type, "replicat_index": rep})
        d = df.loc[df.type.isin([1,2])]
        X = np.array([d["x"], d["y"], d["z"]]).T
        d["sill"] = silhouette_samples(X, labels = d["type"])
        d["chs"] = calinski_harabasz_score(X, labels = d["type"])
        d["cs"] = cs
        # y.append(d.groupby("type",as_index=False).mean())
        y.append(d)

df = pd.concat(y)
#%%
factor = 1.03
rc_ticks['figure.figsize'] = [1.72 * factor, 1.386 * factor]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
sns.lineplot(x = "cs", y = "sill", data=df.loc[df.type == 2], errorbar="se", color="black")
# sns.lineplot(x = "cs", y = "chs", data=df.loc[df.type == 2], errorbar="se", color="black")
# plt.legend(["background","tregs"])
plt.xlim(0,1)
plt.xticks([0, 0.5, 1])
plt.ylim(0, 0.5)
plt.yticks([0, 0.25, 0.5])
plt.xlabel(r"$\varphi$")
plt.ylabel("silhouette score")
plt.tight_layout()
plt.savefig(r"/home/brunner/Documents/Current work/2023_12_08/silhouette_score_over_cs.pdf", bbox_inches='tight', transparent=True)
plt.show()
#
#%%
factor = 1.03
rc_ticks['figure.figsize'] = [1.72 * factor, 1.386 * factor]
sns.set_theme(context="talk", style="ticks", rc=rc_ticks)
fig, ax = plt.subplots()
# sns.lineplot(x = "cs", y = "sill", data=df.loc[df.type == 2], errorbar="se", color="black")
sns.lineplot(x = "cs", y = "chs", data=df.loc[df.type == 2], errorbar="se", color="black")
# plt.legend(["background","tregs"])
plt.xlim(0,1)
plt.xticks([0, 0.5, 1])
plt.ylim(0, 100)
plt.yticks([0, 50, 100])
plt.xlabel(r"$\varphi$")
plt.ylabel("Calinski Harabasz score")
plt.tight_layout()
plt.savefig(r"/home/brunner/Documents/Current work/2023_12_08/calinski_harabasz_score_over_cs.pdf", bbox_inches='tight', transparent=True)
plt.show()
#
