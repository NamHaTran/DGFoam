import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# ============================================================
# Owner and Neighbor cell vertices
# ============================================================

owner_verts = np.array([
    [0.00409653, -0.00341803, 0.0],
    [0.760894,    0.0626306,  0.414656],
    [0.282173,    0.458583,   0.0],
    [0.900843,    0.434146,   0.780304]
])

neighbor_verts = np.array([
    [0.900843, 0.434146, 0.780304],
    [0.760894, 0.0626306, 0.414656],
    [0.282173, 0.458583, 0.0],
    [0.939652, 0.342131, 0.34994]
])


# ============================================================
# Face vertices (CCW from owner)
# ============================================================

P0 = owner_verts[3]
P1 = owner_verts[1]
P2 = owner_verts[2]

face_vertices = np.array([P0, P1, P2])


# ============================================================
# 2x2 Gauss on square
# ============================================================

gauss2D = np.array([
    [-0.57735, -0.57735],
    [ 0.57735, -0.57735],
    [-0.57735,  0.57735],
    [ 0.57735,  0.57735]
])


def square_to_triangle(xi, eta):
    a = 0.5*(1 + xi)
    b = 0.5*(1 + eta)
    r = a
    s = b*(1 - a)
    return r, s


def map_to_physical(r, s, P0, P1, P2):
    return (1 - r - s)*P0 + r*P1 + s*P2


# ============================================================
# Collapsed tet inverse mapping
# ============================================================

def mapXtoEta_tet_collapsed(X, verts, eps=1e-14):

    A, B, C, D = verts

    M = np.column_stack((B - A, C - A, D - A))
    rhs = X - A

    if abs(np.linalg.det(M)) < eps:
        raise ValueError("Degenerate tetra")

    eta_lin = np.linalg.solve(M, rhs)

    lamA = 1 - eta_lin.sum()
    lamB = eta_lin[0]
    lamC = eta_lin[1]
    lamD = eta_lin[2]

    eta3 = 2*lamD - 1
    eta2 = 2*lamC/(1 - lamD) - 1
    eta1 = 2*lamB/(1 - lamC - lamD) - 1

    return np.array([eta1, eta2, eta3])


# ============================================================
# Compute physical Gauss + eta
# ============================================================

physical_gauss = []
eta_owner = []
eta_neighbor = []

for xi, eta in gauss2D:
    r, s = square_to_triangle(xi, eta)
    X = map_to_physical(r, s, P0, P1, P2)

    physical_gauss.append(X)
    eta_owner.append(mapXtoEta_tet_collapsed(X, owner_verts))
    eta_neighbor.append(mapXtoEta_tet_collapsed(X, neighbor_verts))

physical_gauss = np.array(physical_gauss)
eta_owner = np.array(eta_owner)
eta_neighbor = np.array(eta_neighbor)

print("eta_owner:")
print(eta_owner)

print("eta_neighbor:")
print(eta_neighbor)


# ============================================================
# Plot
# ============================================================

fig = plt.figure(figsize=(12,6))

# ---------------- Physical space ----------------
ax1 = fig.add_subplot(121, projection='3d')

tri = Poly3DCollection([face_vertices], alpha=0.4)
ax1.add_collection3d(tri)

ax1.scatter(physical_gauss[:,0],
            physical_gauss[:,1],
            physical_gauss[:,2],
            color='blue', s=60)

ax1.set_title("Physical Face & Gauss Points")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
ax1.set_zlabel("Z")

# ---------------- Reference hex ----------------
ax2 = fig.add_subplot(122, projection='3d')

# Draw reference cube
cube = np.array([
    [-1,-1,-1],[ 1,-1,-1],[ 1, 1,-1],[-1, 1,-1],
    [-1,-1, 1],[ 1,-1, 1],[ 1, 1, 1],[-1, 1, 1]
])

edges = [(0,1),(1,2),(2,3),(3,0),
         (4,5),(5,6),(6,7),(7,4),
         (0,4),(1,5),(2,6),(3,7)]

for e in edges:
    ax2.plot(*zip(cube[e[0]], cube[e[1]]), color='gray')

# Plot eta points
ax2.scatter(eta_owner[:,0],
            eta_owner[:,1],
            eta_owner[:,2],
            color='red', s=60, label='Owner eta')

ax2.scatter(eta_neighbor[:,0],
            eta_neighbor[:,1],
            eta_neighbor[:,2],
            color='black', s=60, label='Neighbor eta')

ax2.set_title("Reference Hex Space")
ax2.set_xlabel("eta1")
ax2.set_ylabel("eta2")
ax2.set_zlabel("eta3")
ax2.legend()

plt.tight_layout()
plt.show()