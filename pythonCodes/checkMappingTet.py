import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.lines import Line2D

# ============================================================
# 1D Gauss points (2-point Gauss-Legendre)
# ============================================================
eta1D = np.array([-0.57735, 0.57735])
n1D = len(eta1D)

# ============================================================
# Generate Gauss points on ref-hex faces (ref coords in [-1,1]^3)
# Face naming:
#   ABCD : xi3 = -1
#   EFGH : xi3 = +1
#   ABEF : xi2 = -1
#   CDGH : xi2 = +1
#   ACEG : xi1 = -1
#   BDFH : xi1 = +1
# ============================================================
faces_hex = {
    "ABCD": [],
    "EFGH": [],
    "ABEF": [],
    "CDGH": [],
    "ACEG": [],
    "BDFH": []
}

for j in range(n1D):
    for i in range(n1D):
        eta1 = float(eta1D[i])
        eta2 = float(eta1D[j])

        faces_hex["ABCD"].append([eta1, eta2, -1.0])
        faces_hex["EFGH"].append([eta1, eta2,  1.0])
        faces_hex["ABEF"].append([eta1, -1.0, eta2])
        faces_hex["CDGH"].append([eta1,  1.0, eta2])
        faces_hex["ACEG"].append([-1.0, eta1, eta2])
        faces_hex["BDFH"].append([ 1.0, eta1, eta2])

# ============================================================
# Ref-hex vertices (A..H) in [-1,1]^3
# ============================================================
V_hex = np.array([
    [-1, -1, -1],  # A 0
    [ 1, -1, -1],  # B 1
    [ 1,  1, -1],  # C 2
    [-1,  1, -1],  # D 3
    [-1, -1,  1],  # E 4
    [ 1, -1,  1],  # F 5
    [ 1,  1,  1],  # G 6
    [-1,  1,  1],  # H 7
], dtype=float)
labels_hex = ["A","B","C","D","E","F","G","H"]

edges_hex = [
    (0,1),(1,2),(2,3),(3,0),  # bottom
    (4,5),(5,6),(6,7),(7,4),  # top
    (0,4),(1,5),(2,6),(3,7)   # verticals
]

faces_hex_polys = {
    "ABCD": [0,1,2,3],
    "EFGH": [4,5,6,7],
    "ABEF": [0,1,5,4],
    "CDGH": [3,2,6,7],
    "ACEG": [0,3,7,4],
    "BDFH": [1,2,6,5],
}

# ============================================================
# Physical tetra vertices (given order -> A,B,C,D)
# ============================================================
V_tet = np.array([
    [0.900843, 0.434146, 0.780304],   # A
    [0.760894, 0.0626306, 0.414656],  # B
    [0.282173, 0.458583, 0.0],        # C
    [0.939652, 0.342131, 0.34994],    # D
], dtype=float)
labels_tet = ["A","B","C","D"]

faces_tet_tris = [
    [0,1,2],  # ABC
    [0,1,3],  # ABD
    [0,2,3],  # ACD
    [1,2,3],  # BCD
]
edges_tet = [(0,1),(1,2),(2,0),(0,3),(1,3),(2,3)]

# ============================================================
# Mapping using the tetra "collapsed-cube" shape functions written in eta
#
# eta = (eta1, eta2, eta3) in [-1,1]^3
#
#   NA(eta) = (1-eta1)/2 * (1-eta2)/2 * (1-eta3)/2
#   NB(eta) = (1+eta1)/2 * (1-eta2)/2 * (1-eta3)/2
#   NC(eta) = (1+eta2)/2 * (1-eta3)/2
#   ND(eta) = (1+eta3)/2
#
# X = NA*A + NB*B + NC*C + ND*D
# ============================================================
def hex_to_tet_physical_shape_eta(xi, Vt):
    eta1, eta2, eta3 = np.asarray(xi, dtype=float)

    NA = 0.5*(1.0 - eta1) * 0.5*(1.0 - eta2) * 0.5*(1.0 - eta3)
    NB = 0.5*(1.0 + eta1) * 0.5*(1.0 - eta2) * 0.5*(1.0 - eta3)
    NC = 0.5*(1.0 + eta2) * 0.5*(1.0 - eta3)
    ND = 0.5*(1.0 + eta3)

    # Sanity check: must sum to 1
    # s = NA + NB + NC + ND
    # if abs(s - 1.0) > 1e-12:
    #     raise RuntimeError(f"Shape sum != 1: {s}")

    return NA*Vt[0] + NB*Vt[1] + NC*Vt[2] + ND*Vt[3]

mapped_tet = {
    k: np.array([hex_to_tet_physical_shape_eta(p, V_tet) for p in pts])
    for k, pts in faces_hex.items()
}

# ============================================================
# Plot helpers
# ============================================================
def set_axes_equal(ax):
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    x_mid = np.mean(x_limits)
    y_mid = np.mean(y_limits)
    z_mid = np.mean(z_limits)

    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_mid - plot_radius, x_mid + plot_radius])
    ax.set_ylim3d([y_mid - plot_radius, y_mid + plot_radius])
    ax.set_zlim3d([z_mid - plot_radius, z_mid + plot_radius])

# Consistent colors per face group
face_order = ["ABCD","EFGH","ABEF","CDGH","ACEG","BDFH"]
cmap = plt.get_cmap("tab10")
colors = {name: cmap(i) for i, name in enumerate(face_order)}

legend_handles = [
    Line2D([0],[0], marker='o', linestyle='None', markersize=7, color=colors[name], label=name)
    for name in face_order
]

# ============================================================
# Figure 1: Reference hex + face Gauss points
# ============================================================
fig1 = plt.figure(figsize=(9, 7))
ax1 = fig1.add_subplot(111, projection="3d")
ax1.set_title("Reference hex ([-1,1]^3) with face Gauss points")

for i, j in edges_hex:
    p0, p1 = V_hex[i], V_hex[j]
    ax1.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]])

for fname, vidx in faces_hex_polys.items():
    poly = [V_hex[idx] for idx in vidx]
    ax1.add_collection3d(Poly3DCollection([poly], alpha=0.07))

for idx, (x,y,z) in enumerate(V_hex):
    ax1.text(x, y, z, f" {labels_hex[idx]}", fontsize=10)

for fname in face_order:
    P = np.array(faces_hex[fname], dtype=float)
    ax1.scatter(P[:,0], P[:,1], P[:,2], s=50, color=colors[fname])

ax1.set_xlabel(r"$\xi_1$")
ax1.set_ylabel(r"$\xi_2$")
ax1.set_zlabel(r"$\xi_3$")
ax1.legend(handles=legend_handles, loc="upper left")
ax1.view_init(elev=20, azim=35)
set_axes_equal(ax1)

# ============================================================
# Figure 2: Physical tetra + mapped Gauss points
# ============================================================
fig2 = plt.figure(figsize=(9, 7))
ax2 = fig2.add_subplot(111, projection="3d")
ax2.set_title("Physical tetra with mapped Gauss points")

tet_polys = [[V_tet[i] for i in tri] for tri in faces_tet_tris]
ax2.add_collection3d(Poly3DCollection(tet_polys, alpha=0.12))

for i, j in edges_tet:
    p0, p1 = V_tet[i], V_tet[j]
    ax2.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]])

for idx, (x,y,z) in enumerate(V_tet):
    ax2.text(x, y, z, f" {labels_tet[idx]}", fontsize=11)

for fname in face_order:
    P = mapped_tet[fname]
    ax2.scatter(P[:,0], P[:,1], P[:,2], s=55, color=colors[fname])

ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("z")
ax2.legend(handles=legend_handles, loc="upper left")
ax2.view_init(elev=20, azim=35)
set_axes_equal(ax2)

plt.show()
