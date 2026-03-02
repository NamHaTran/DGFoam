import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# ============================================================
# 1D Gauss points (2-point Gauss-Legendre)
# ============================================================
eta1D = np.array([-0.57735, 0.57735])
n1D = len(eta1D)


# ============================================================
# Generate Gauss points on ref-hex faces (ref coords in [-1,1]^3)
# Face naming follows your convention:
#   ABCD : xi3 = -1
#   EFGH : xi3 = +1
#   ABEF : xi2 = -1
#   CDGH : xi2 = +1
#   ACEG : xi1 = -1
#   BDFH : xi1 = +1
# ============================================================
faces_ref = {
    "ABCD": [],
    "EFGH": [],
    "ABEF": [],
    "CDGH": [],
    "ACEG": [],
    "BDFH": []
}

for j in range(n1D):
    for i in range(n1D):
        xi1 = eta1D[i]
        xi2 = eta1D[j]

        faces_ref["ABCD"].append([xi1, xi2, -1.0])
        faces_ref["EFGH"].append([xi1, xi2,  1.0])
        faces_ref["ABEF"].append([xi1, -1.0, xi2])
        faces_ref["CDGH"].append([xi1,  1.0, xi2])
        faces_ref["ACEG"].append([-1.0, xi1, xi2])
        faces_ref["BDFH"].append([ 1.0, xi1, xi2])

for k in faces_ref:
    faces_ref[k] = np.array(faces_ref[k], dtype=float)


# ============================================================
# Ref hex vertices (standard) in order [A,B,C,D,E,F,G,H]
# ============================================================
Vref = np.array([
    [-1.0, -1.0, -1.0],  # A (0)
    [ 1.0, -1.0, -1.0],  # B (1)
    [ 1.0,  1.0, -1.0],  # C (2)
    [-1.0,  1.0, -1.0],  # D (3)
    [-1.0, -1.0,  1.0],  # E (4)
    [ 1.0, -1.0,  1.0],  # F (5)
    [ 1.0,  1.0,  1.0],  # G (6)
    [-1.0,  1.0,  1.0],  # H (7)
], dtype=float)

vNames = ["A", "B", "C", "D", "E", "F", "G", "H"]


# ============================================================
# Physical hex vertices (your data) in the SAME order [A..H]
# ((0.05 0 0.25) (0.05 0 0.5) (0.05 0.25 0.5) (0.05 0.25 0.25)
#  (0 0 0.25) (0 0 0.5) (0 0.25 0.5) (0 0.25 0.25))
# ============================================================
Vphys = np.array([
    [0.05, 0.00, 0.25],  # A
    [0.05, 0.00, 0.50],  # B
    [0.05, 0.25, 0.50],  # C
    [0.05, 0.25, 0.25],  # D
    [0.00, 0.00, 0.25],  # E
    [0.00, 0.00, 0.50],  # F
    [0.00, 0.25, 0.50],  # G
    [0.00, 0.25, 0.25],  # H
], dtype=float)


# ============================================================
# Hex trilinear shape functions with requested swaps:
#   - swap shape basis of vertex C <-> D
#   - swap shape basis of vertex G <-> H
#
# Returned order: [A,B,C,D,E,F,G,H]
# ============================================================
def hex_shape_swapped(xi):
    xi1, xi2, xi3 = xi

    a1m = 0.5 * (1.0 - xi1)
    a1p = 0.5 * (1.0 + xi1)
    a2m = 0.5 * (1.0 - xi2)
    a2p = 0.5 * (1.0 + xi2)
    a3m = 0.5 * (1.0 - xi3)
    a3p = 0.5 * (1.0 + xi3)

    # Standard shapes
    NA = a1m * a2m * a3m
    NB = a1p * a2m * a3m
    NC = a1p * a2p * a3m 
    ND = a1m * a2p * a3m

    NE = a1m * a2m * a3p
    NF = a1p * a2m * a3p
    NG = a1p * a2p * a3p
    NH = a1m * a2p * a3p

    return np.array([NA, NB, NC, ND, NE, NF, NG, NH], dtype=float)


def map_refhex_to_physhex(xi, Xphys):
    N = hex_shape_swapped(xi)  # (8,)
    return (N[:, None] * Xphys).sum(axis=0)


# ============================================================
# Map all Gauss points from ref hex -> physical hex
# ============================================================
faces_phys = {}
for fname, pts in faces_ref.items():
    faces_phys[fname] = np.array([map_refhex_to_physhex(p, Vphys) for p in pts], dtype=float)


# ============================================================
# Geometry topology (for drawing)
# ============================================================
hexFaces_vid = [
    [0, 1, 2, 3],  # ABCD
    [4, 5, 6, 7],  # EFGH
    [0, 1, 5, 4],  # ABEF
    [3, 2, 6, 7],  # DCGH (same quad as CDGH)
    [0, 3, 7, 4],  # ADHE (left side quad)
    [1, 2, 6, 5],  # BCFG (right side quad)
]

hexEdges_vid = [
    (0, 1), (1, 2), (2, 3), (3, 0),
    (4, 5), (5, 6), (6, 7), (7, 4),
    (0, 4), (1, 5), (2, 6), (3, 7)
]

gaussColors = {
    "ABCD": "red",
    "EFGH": "blue",
    "ABEF": "green",
    "CDGH": "cyan",
    "ACEG": "black",
    "BDFH": "yellow"
}


# ============================================================
# Plot helpers
# ============================================================
def add_hex_geometry_wire(ax, V, title):
    # Faces (transparent)
    facePolys = [[V[i] for i in f] for f in hexFaces_vid]
    poly = Poly3DCollection(facePolys, alpha=0.12, edgecolor="k", linewidths=0.8)
    ax.add_collection3d(poly)

    # Edges
    for (i, j) in hexEdges_vid:
        p = V[i]
        q = V[j]
        ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], "k", linewidth=1.2)

    # Vertices + labels
    ax.scatter(V[:, 0], V[:, 1], V[:, 2], s=60, color="black")
    for k, name in enumerate(vNames):
        ax.text(V[k, 0], V[k, 1], V[k, 2], f"  {name}", fontsize=11)

    # Axis and title
    ax.set_title(title)
    ax.set_box_aspect([1, 1, 1])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")


def add_hex_geometry_filled(ax, V, title):
    # Faces (filled)
    facePolys = [[V[i] for i in f] for f in hexFaces_vid]
    poly = Poly3DCollection(
        facePolys,
        facecolor="lightsteelblue",
        edgecolor="k",
        linewidths=1.0,
        alpha=0.35
    )
    ax.add_collection3d(poly)

    # Edges
    for (i, j) in hexEdges_vid:
        p = V[i]
        q = V[j]
        ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], "k", linewidth=1.5)

    # Vertices + labels
    ax.scatter(V[:, 0], V[:, 1], V[:, 2], s=60, color="black")
    for k, name in enumerate(vNames):
        ax.text(V[k, 0], V[k, 1], V[k, 2], f"  {name}", fontsize=11)

    # Fix axis limits to cell bounding box
    xmin, ymin, zmin = V.min(axis=0)
    xmax, ymax, zmax = V.max(axis=0)

    pad = 0.02 * max(xmax - xmin, ymax - ymin, zmax - zmin)
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(ymin - pad, ymax + pad)
    ax.set_zlim(zmin - pad, zmax + pad)

    ax.set_box_aspect([xmax - xmin, ymax - ymin, zmax - zmin])

    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")


def add_gauss_points(ax, facesPtsDict):
    for fname, pts in facesPtsDict.items():
        ax.scatter(
            pts[:, 0], pts[:, 1], pts[:, 2],
            color=gaussColors[fname],
            s=70,
            label=fname
        )


# ============================================================
# Figure 1: Ref hex + Gauss points on ref faces
# ============================================================
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection="3d")
add_hex_geometry_wire(ax1, Vref, "Ref Hex ([-1,1]^3) with Face Gauss Points")
add_gauss_points(ax1, faces_ref)
ax1.legend()


# ============================================================
# Figure 2: Physical hex (filled) + mapped Gauss points
# ============================================================
fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection="3d")
add_hex_geometry_filled(ax2, Vphys, "Physical Hex (Filled) with Mapped Face Gauss Points (Swapped Shapes)")
add_gauss_points(ax2, faces_phys)
ax2.legend()

plt.show()