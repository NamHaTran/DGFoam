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
# Generate Gauss points on ref hex faces
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

        eta1 = eta1D[i]
        eta2 = eta1D[j]

        faces_hex["ABCD"].append([eta1, eta2, -1.0])
        faces_hex["EFGH"].append([eta1, eta2,  1.0])
        faces_hex["ABEF"].append([eta1, -1.0, eta2])
        faces_hex["CDGH"].append([eta1,  1.0, eta2])
        faces_hex["ACEG"].append([-1.0, eta1, eta2])
        faces_hex["BDFH"].append([ 1.0, eta1, eta2])

for k in faces_hex:
    faces_hex[k] = np.array(faces_hex[k], dtype=float)

# ============================================================
# Ref hex vertices (A..H) in [-1,1]^3
# ============================================================
Vhex = {
    "A": np.array([-1.0, -1.0, -1.0]),
    "B": np.array([ 1.0, -1.0, -1.0]),
    "C": np.array([ 1.0,  1.0, -1.0]),
    "D": np.array([-1.0,  1.0, -1.0]),
    "E": np.array([-1.0, -1.0,  1.0]),
    "F": np.array([ 1.0, -1.0,  1.0]),
    "G": np.array([ 1.0,  1.0,  1.0]),
    "H": np.array([-1.0,  1.0,  1.0]),
}

# Hex edges (for plotting)
hex_edges = [
    ("A","B"), ("B","C"), ("C","D"), ("D","A"),
    ("E","F"), ("F","G"), ("G","H"), ("H","E"),
    ("A","E"), ("B","F"), ("C","G"), ("D","H")
]

# Hex faces as vertex loops
hex_face_verts = {
    "ABCD": ["A","B","C","D"],
    "EFGH": ["E","F","G","H"],
    "ABEF": ["A","B","F","E"],
    "CDGH": ["C","D","H","G"],
    "ACEG": ["A","D","H","E"],  # x=-1 face (A-D-H-E)
    "BDFH": ["B","C","G","F"],  # x=+1 face (B-C-G-F)
}

# ============================================================
# Physical prism vertices (6 points) - label as A..F
# (Bottom triangle: A,B,C ; Top triangle: D,E,F)
# ============================================================
Vpr = {
    "A": np.array([17.5629, 72.7914, 23.2775]),
    "B": np.array([18.4082, 75.4340, 23.8272]),
    "C": np.array([19.4915, 74.8623, 22.2462]),
    "D": np.array([16.2368, 73.4544, 22.1290]),
    "E": np.array([17.0820, 76.0970, 22.6787]),
    "F": np.array([18.1654, 75.5254, 21.0977]),
}

# Prism edges (wedge): bottom tri + top tri + verticals
pr_edges = [
    ("A","B"), ("B","C"), ("C","A"),
    ("D","E"), ("E","F"), ("F","D"),
    ("A","D"), ("B","E"), ("C","F")
]

# Prism faces (triangles and quads)
pr_faces = {
    "ABC": ["A","B","C"],     # bottom tri
    "DEF": ["D","E","F"],     # top tri
    "ABED": ["A","B","E","D"],# quad
    "BCFE": ["B","C","F","E"],# quad
    "CADF": ["C","A","D","F"] # quad
}

# ============================================================
# Mapping: ref-hex point (xi1,xi2,xi3) -> physical wedge (A..F)
# Use Duffy (square->triangle) for (xi1,xi2), linear for xi3.
# ============================================================
def duffy_square_to_triangle(u, v):
    """
    Map (u,v) in [-1,1]^2 (square) to (a,b) in triangle coordinates (a>=0, b>=0, a+b<=1).
    a corresponds to "towards vertex B", b corresponds to "towards vertex C".
    """
    # map to [0,1]
    uu = 0.5*(u + 1.0)
    vv = 0.5*(v + 1.0)

    a = uu
    b = vv*(1.0 - a)
    return a, b

def wedge_map_from_hex_ref(xi):
    """
    Map xi=(xi1,xi2,xi3) from ref hex [-1,1]^3 to physical prism using:
      - (xi1,xi2) -> triangle barycentric weights via Duffy
      - xi3 -> interpolation between bottom and top triangles
    """
    xi1, xi2, xi3 = xi

    a, b = duffy_square_to_triangle(xi1, xi2)
    N1 = 1.0 - a - b
    N2 = a
    N3 = b

    # bottom triangle (A,B,C) and top triangle (D,E,F)
    Pb = N1*Vpr["A"] + N2*Vpr["B"] + N3*Vpr["C"]
    Pt = N1*Vpr["D"] + N2*Vpr["E"] + N3*Vpr["F"]

    lam = 0.5*(xi3 + 1.0)  # xi3=-1 -> 0 (bottom), xi3=+1 -> 1 (top)
    X = (1.0 - lam)*Pb + lam*Pt
    return X

# Map all face Gauss points
faces_prism = {name: np.array([wedge_map_from_hex_ref(p) for p in pts]) for name, pts in faces_hex.items()}

# ============================================================
# Plot helpers
# ============================================================
def add_edges(ax, Vdict, edges, lw=1.5, alpha=1.0):
    for a, b in edges:
        P = Vdict[a]
        Q = Vdict[b]
        ax.plot([P[0], Q[0]], [P[1], Q[1]], [P[2], Q[2]], linewidth=lw, alpha=alpha)

def add_vertex_labels(ax, Vdict, offset=(0.0, 0.0, 0.0), fontsize=10):
    dx, dy, dz = offset
    for name, p in Vdict.items():
        ax.text(p[0]+dx, p[1]+dy, p[2]+dz, name, fontsize=fontsize)

def add_faces(ax, Vdict, face_loops, alpha=0.10):
    polys = []
    for _, loop in face_loops.items():
        polys.append([Vdict[v] for v in loop])
    coll = Poly3DCollection(polys, alpha=alpha)
    ax.add_collection3d(coll)

def equal_aspect_3d(ax, pts):
    pts = np.asarray(pts)
    mins = pts.min(axis=0)
    maxs = pts.max(axis=0)
    ctr = 0.5*(mins + maxs)
    span = (maxs - mins).max()
    if span <= 0:
        span = 1.0
    r = 0.5*span
    ax.set_xlim(ctr[0]-r, ctr[0]+r)
    ax.set_ylim(ctr[1]-r, ctr[1]+r)
    ax.set_zlim(ctr[2]-r, ctr[2]+r)

# Face colors (consistent across both figures)
face_colors = {
    "ABCD": "tab:blue",
    "EFGH": "tab:orange",
    "ABEF": "tab:green",
    "CDGH": "tab:red",
    "ACEG": "tab:purple",
    "BDFH": "tab:brown",
}

# ============================================================
# Figure 1: Reference hex + all face Gauss points
# ============================================================
fig1 = plt.figure(figsize=(10, 7))
ax1 = fig1.add_subplot(111, projection="3d")
ax1.set_title("Reference Hex [-1,1]^3: vertices/edges/faces + Gauss points on faces")

# Draw hex faces (transparent) and edges
add_faces(ax1, Vhex, hex_face_verts, alpha=0.08)
add_edges(ax1, Vhex, hex_edges, lw=1.8, alpha=0.9)

# Labels
add_vertex_labels(ax1, Vhex, offset=(0.03, 0.03, 0.03), fontsize=11)

# Plot Gauss points per face with color
for fname, pts in faces_hex.items():
    c = face_colors.get(fname, "k")
    ax1.scatter(pts[:,0], pts[:,1], pts[:,2], s=55, color=c, depthshade=True)

# Legend
legend_handles = [Line2D([0],[0], marker='o', color='w', label=f, markerfacecolor=face_colors[f], markersize=10)
                  for f in faces_hex.keys()]
ax1.legend(handles=legend_handles, loc="upper left")

# Axes & aspect
all_hex_pts = np.array(list(Vhex.values()))
equal_aspect_3d(ax1, all_hex_pts)
ax1.set_xlabel("xi1")
ax1.set_ylabel("xi2")
ax1.set_zlabel("xi3")

# ============================================================
# Figure 2: Physical prism + mapped Gauss points (same colors)
# ============================================================
fig2 = plt.figure(figsize=(10, 7))
ax2 = fig2.add_subplot(111, projection="3d")
ax2.set_title("Physical Prism: vertices/edges/faces + mapped Gauss points (same colors as ref hex)")

# Draw prism faces (transparent) and edges
add_faces(ax2, Vpr, pr_faces, alpha=0.10)
add_edges(ax2, Vpr, pr_edges, lw=1.8, alpha=0.9)

# Labels
# Slight offset to reduce overlap with markers
add_vertex_labels(ax2, Vpr, offset=(0.03, 0.03, 0.03), fontsize=11)

# Plot mapped Gauss points per original hex face
for fname, pts in faces_prism.items():
    c = face_colors.get(fname, "k")
    ax2.scatter(pts[:,0], pts[:,1], pts[:,2], s=55, color=c, depthshade=True)

# Legend (same)
ax2.legend(handles=legend_handles, loc="upper left")

# Axes & aspect
all_pr_pts = np.array(list(Vpr.values()))
equal_aspect_3d(ax2, all_pr_pts)
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("z")

plt.show()
