import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# 1. Trilinear shape functions for HEX
#
# Reference cube: xi1, xi2, xi3 in [-1,1]
#
# Node ordering (OpenFOAM-style, consistent & important):
#
#   bottom (xi3 = -1):
#     0: (-1,-1,-1)  A
#     1: ( 1,-1,-1)  B
#     2: ( 1, 1,-1)  C
#     3: (-1, 1,-1)  D
#
#   top (xi3 = +1):
#     4: (-1,-1, 1)  A'
#     5: ( 1,-1, 1)  B'
#     6: ( 1, 1, 1)  C'
#     7: (-1, 1, 1)  D'
# ============================================================

def hex_shape_functions(xi):
    xi1, xi2, xi3 = xi

    N = np.zeros(8)

    signs = [
        (-1,-1,-1),
        ( 1,-1,-1),
        ( 1, 1,-1),
        (-1, 1,-1),
        (-1,-1, 1),
        ( 1,-1, 1),
        ( 1, 1, 1),
        (-1, 1, 1),
    ]

    for a, (i, j, k) in enumerate(signs):
        N[a] = 0.125 * (1 + i*xi1) * (1 + j*xi2) * (1 + k*xi3)

    return N


# ============================================================
# 2. Explicit mapping x(xi)
# ============================================================

def x_of_xi_hex(xi, X8):
    """
    Explicit mapping x(xi) for straight-sided HEX.

    X8: array (8,3) of physical hex vertices
    """
    N = hex_shape_functions(xi)
    x = (N[:, None] * X8).sum(axis=0)
    return x


# ============================================================
# 3. Plot helpers
# ============================================================

def hex_edges():
    return [
        (0,1),(1,2),(2,3),(3,0),   # bottom
        (4,5),(5,6),(6,7),(7,4),   # top
        (0,4),(1,5),(2,6),(3,7)    # vertical
    ]


def plot_hex_wire(ax, P):
    for i, j in hex_edges():
        ax.plot([P[i,0], P[j,0]],
                [P[i,1], P[j,1]],
                [P[i,2], P[j,2]])


# ============================================================
# 4. Main
# ============================================================

def main():

    # --------------------------------------------------------
    # Physical HEX vertices
    # --------------------------------------------------------
    X8 = np.array([
        [0.0, 0.0, 0.0],   # A
        [2.0, 0.2, 0.0],   # B
        [2.2, 1.8, 0.1],   # C
        [0.1, 1.7, 0.0],   # D
        [0.0, 0.0, 1.5],   # A'
        [2.1, 0.3, 1.6],   # B'
        [2.3, 1.9, 1.7],   # C'
        [0.2, 1.8, 1.6],   # D'
    ])

    # --------------------------------------------------------
    # Interior points in reference cube
    # --------------------------------------------------------
    npts = 300
    xi_pts = np.random.uniform(-1.0, 1.0, size=(npts, 3))

    x_phys = np.zeros((npts, 3))
    for i in range(npts):
        x_phys[i] = x_of_xi_hex(xi_pts[i], X8)

    # --------------------------------------------------------
    # Plot physical HEX
    # --------------------------------------------------------
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection="3d")
    plot_hex_wire(ax1, X8)
    ax1.scatter(X8[:,0], X8[:,1], X8[:,2], s=60, label="Nodes")
    ax1.scatter(x_phys[:,0], x_phys[:,1], x_phys[:,2],
                s=10, label="x(xi)")
    ax1.set_title("Physical HEX: explicit x(xi)")
    ax1.legend()

    # --------------------------------------------------------
    # Plot reference cube
    # --------------------------------------------------------
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection="3d")
    ax2.scatter(xi_pts[:,0], xi_pts[:,1], xi_pts[:,2], s=10)
    ax2.set_xlim([-1,1])
    ax2.set_ylim([-1,1])
    ax2.set_zlim([-1,1])
    ax2.set_title("Reference cube (xi)")

    plt.show()


if __name__ == "__main__":
    main()
