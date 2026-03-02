import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ============================================================
# Physical tetra vertices (Owner cell)
# ============================================================
A = np.array([0.00409653, -0.00341803, 0.0])
B = np.array([0.760894,    0.0626306,  0.414656])
C = np.array([0.282173,    0.458583,   0.0])
D = np.array([0.900843,    0.434146,   0.780304])

vertices_phys = {"A": A, "B": B, "C": C, "D": D}

# ============================================================
# 1D Gauss points
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

# ============================================================
# hex -> ref tetra mapping (Duffy)
# ============================================================
def hex_to_tet(eta):
    eta1, eta2, eta3 = eta

    xi3 = eta3
    xi2 = 0.5 * (1 - eta3) * (1 + eta2) - 1.0
    xi1 = 0.5 * (1 - eta3) * (1 + eta1) - 1.0

    return np.array([xi1, xi2, xi3])

# ============================================================
# ref tetra -> physical tetra (shape functions)
# ============================================================
def ref_to_phys_tet(eta):

    e1, e2, e3 = eta

    NA = 0.5*(1-e1) * 0.5*(1-e2) * 0.5*(1-e3)
    NB = 0.5*(1+e1) * 0.5*(1-e2) * 0.5*(1-e3)
    NC = 0.5*(1-e1) * 0.5*(1+e2) * 0.5*(1-e3)
    ND = 0.5*(1+e3)

    return NA*A + NB*B + NC*C + ND*D

# ============================================================
# Map all Gauss points to physical tetra
# ============================================================
faces_phys = {}

for name, pts in faces_hex.items():

    mapped = []
    for p in pts:
        #eta_tet = hex_to_tet(p)
        phys_pt = ref_to_phys_tet(p)
        mapped.append(phys_pt)

    faces_phys[name] = np.array(mapped)

# ============================================================
# Plot physical tetra
# ============================================================
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot edges
edges = [(A,B),(A,C),(A,D),(B,C),(B,D),(C,D)]
for p,q in edges:
    ax.plot([p[0],q[0]],
            [p[1],q[1]],
            [p[2],q[2]], 'k')

# Plot vertices + labels
for name, coord in vertices_phys.items():
    ax.scatter(coord[0], coord[1], coord[2],
               color='magenta', s=120)
    ax.text(coord[0], coord[1], coord[2],
            f"  {name}",
            fontsize=12, color='magenta')

# ============================================================
# Plot Gauss points (same color rule)
# ============================================================
colors = {
    "ABCD": "red",
    "EFGH": "blue",
    "ABEF": "green",
    "CDGH": "cyan",
    "ACEG": "black",
    "BDFH": "yellow"
}

for name, pts in faces_phys.items():
    ax.scatter(pts[:,0], pts[:,1], pts[:,2],
               color=colors[name], label=name, s=60)

ax.set_box_aspect([1,1,1])
ax.legend()
plt.title("Physical Tetra with Mapped Gauss Points")
plt.show()