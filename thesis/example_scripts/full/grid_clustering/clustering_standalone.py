import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from thesis.cellBehaviourUtilities.grid_clustering import get_cell_grid_positions, get_apc_positions, make_clusters
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from thesis.cellBehaviourUtilities.grid_clustering import get_cell_grid_positions, get_apc_positions, make_clusters

"""
This script demonstrates the usage of the grid clustering utilities, without an actual simulation.
 
"""

"""get cell grid. In simulation this would be generated from entity positions."""
cell_grid_positions = get_cell_grid_positions(200, 200, 200)

"""Creates some apcs"""
apcs = get_apc_positions(cell_grid_positions, no_apcs=4)
# apcs = np.array([[50,50,50]])

"""fractions for each cell type and corresponding cluster strengths"""
fractions = [0.05, 0.1]
cluster_strengths = [0, 0]

cell_type = make_clusters(cell_grid_positions, apcs, fractions, cluster_strengths)

df = pd.DataFrame(
    {"x": cell_grid_positions[:, 0], "y": cell_grid_positions[:, 1], "z": cell_grid_positions[:, 2], "type": cell_type})

"""creates 3D plot"""
names = ["apc", "naive", "sec", "treg"]
colors = ["gray", "red", "blue"]
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection='3d')
ax.set_xlim([df["x"].min(), df["x"].max()])
ax.set_ylim([df["y"].min(), df["y"].max()])
ax.set_zlim([df["z"].min(), df["z"].max()])
ax.scatter(apcs[:, 0], apcs[:, 1], apcs[:, 2], s=200, color="green")
for i, d in df.groupby(["type"]):
    if np.isin(i, [0, 1, 2]):
        ax.scatter(d["x"], d["y"], d["z"], s=50, color=colors[int(d["type"].unique())])

plt.legend(names)
plt.show()
print(df.groupby("type").count()["x"] / len(df))
