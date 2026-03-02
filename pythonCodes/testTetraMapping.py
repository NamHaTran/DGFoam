# ---------------------------------------------------------------
# Improved tetra visualization
# - Merge overlapping vertex labels (e.g., Eo/An) into a single label
# - Add legend for: Owner faces/verts/Gauss, Neighbor faces/verts/Gauss
# - Keep colors:
#     Owner faces = cyan
#     Neighbor faces = orange
#     Owner vertices+Gauss = black
#     Neighbor vertices+Gauss = red
# - Equal axis scaling
# 彩森一花
# ---------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# -------------------------
# Data (as provided)
# -------------------------
owner_vertices = np.array([
    [0.00409653, -0.00341803, 0.0],       # Ao
    [0.760894,    0.0626306,  0.414656],  # Bo
    [0.282173,    0.458583,   0.0],       # Co
    [0.900843,    0.434146,   0.780304]   # Eo
])

neighbor_vertices = np.array([
    [0.900843,    0.434146,   0.780304],  # An  (same as Eo)
    [0.760894,    0.0626306,  0.414656],  # Bn  (same as Bo)
    [0.282173,    0.458583,   0.0],       # Cn  (same as Co)
    [0.939652,    0.342131,   0.34994]    # En
])

owner_labels = ["Ao", "Bo", "Co", "Eo"]
neighbor_labels = ["An", "Bn", "Cn", "En"]

# Mapping function
def map_eta_to_x_tet(eta, vertices):
    eta1, eta2, eta3 = eta
    NA = 0.125 * (1 - eta1) * (1 - eta2) * (1 - eta3)
    NB = 0.125 * (1 + eta1) * (1 - eta2) * (1 - eta3)
    NC = 0.25 * (1 + eta2) * (1 - eta3)
    NE = 0.5   * (1 + eta3)
    return NA*vertices[0] + NB*vertices[1] + NC*vertices[2] + NE*vertices[3]

# Gauss points (reference) - same as before
owner_eta = np.array([
    [1,         -0.11814611,  0.24401651],
    [1,         -0.89282018, -0.66666651],
    [1,          0.49281995, -0.66666651],
    [1,         -0.65108456, -0.91068349]
])

neighbor_eta = np.array([
    [-0.49281995, -0.66666651, -1],
    [0.65108456, -0.91068349, -1],
    [0.11814611,  0.24401651, -1],
    [0.89282018, -0.66666651, -1]
])

owner_phys = np.array([map_eta_to_x_tet(e, owner_vertices) for e in owner_eta])
neighbor_phys = np.array([map_eta_to_x_tet(e, neighbor_vertices) for e in neighbor_eta])

faces = [
    [0,1,2],
    [0,1,3],
    [0,2,3],
    [1,2,3]
]

# -------------------------
# Build merged labels for coincident vertices
# -------------------------
tol = 1e-10

# Collect all points with their labels
all_pts = []
for p, lab in zip(owner_vertices, owner_labels):
    all_pts.append(("owner", p, lab))
for p, lab in zip(neighbor_vertices, neighbor_labels):
    all_pts.append(("neighbor", p, lab))

# Cluster points by proximity
clusters = []  # list of dict: {"p": representative, "labels": [..]}
for src, p, lab in all_pts:
    placed = False
    for c in clusters:
        if np.linalg.norm(p - c["p"]) < tol:
            c["labels"].append(lab)
            placed = True
            break
    if not placed:
        clusters.append({"p": p.copy(), "labels": [lab]})

# -------------------------
# Plot
# -------------------------
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Faces
for f in faces:
    poly = Poly3DCollection([owner_vertices[f]], alpha=0.3, facecolor='cyan', edgecolor='c')
    ax.add_collection3d(poly)

for f in faces:
    poly = Poly3DCollection([neighbor_vertices[f]], alpha=0.3, facecolor='orange', edgecolor='darkorange')
    ax.add_collection3d(poly)

# Vertices
ax.scatter(owner_vertices[:,0], owner_vertices[:,1], owner_vertices[:,2], color='black', s=60)
ax.scatter(neighbor_vertices[:,0], neighbor_vertices[:,1], neighbor_vertices[:,2], color='red', s=60)

# Gauss points
ax.scatter(owner_phys[:,0], owner_phys[:,1], owner_phys[:,2], color='black', s=40)
ax.scatter(neighbor_phys[:,0], neighbor_phys[:,1], neighbor_phys[:,2], color='red', s=40)

# Label offsets (small relative offset based on plot scale)
plot_pts = np.vstack((owner_vertices, neighbor_vertices, owner_phys, neighbor_phys))
mins = plot_pts.min(axis=0)
maxs = plot_pts.max(axis=0)
max_range = (maxs - mins).max()
centers = (maxs + mins) / 2.0

# Offset magnitude
d = 0.03 * max_range

# Put merged labels once per cluster
# Alternate offset directions so multiple labels don't stack
dirs = np.array([
    [ 1,  1,  1],
    [ 1, -1,  1],
    [-1,  1,  1],
    [-1, -1,  1],
    [ 1,  1, -1],
    [ 1, -1, -1],
    [-1,  1, -1],
    [-1, -1, -1],
], dtype=float)

for k, c in enumerate(clusters):
    p = c["p"]
    labs = "/".join(c["labels"])
    off = d * dirs[k % len(dirs)]
    ax.text(p[0] + off[0], p[1] + off[1], p[2] + off[2], labs, color='k')

# Equal axis scaling
ax.set_xlim(centers[0]-max_range/2, centers[0]+max_range/2)
ax.set_ylim(centers[1]-max_range/2, centers[1]+max_range/2)
ax.set_zlim(centers[2]-max_range/2, centers[2]+max_range/2)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

# Legend (clear ownership)
legend_handles = [
    Patch(facecolor='cyan', edgecolor='c', alpha=0.3, label='Owner cell faces'),
    Patch(facecolor='orange', edgecolor='darkorange', alpha=0.3, label='Neighbor cell faces'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='black', markersize=8, label='Owner vertices/Gauss'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='red', markersize=8, label='Neighbor vertices/Gauss'),
]
ax.legend(handles=legend_handles, loc='upper left')

plt.show()
