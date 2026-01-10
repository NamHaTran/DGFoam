import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# 1. Random points INSIDE reference prism (xi-coordinates)
#    Prism is embedded in reference cube [-1,1]^3
# ============================================================

def random_points_reference_prism(npts):
    """
    Generate random points such that:
      - (xi1, xi3) lie inside the red triangle
      - xi2 in [-1, 1]

    Triangle in (xi1, xi3):
      -1 <= xi3 <= 1
      -1 <= xi1 <= -1 + (1 - xi3)/2
    """
    xi = np.zeros((npts, 3))

    for i in range(npts):
        xi1 = np.random.uniform(-1.0, 1.0)
        xi3 = np.random.uniform(-1.0, -xi1)  # upper bound depends on xi1
        xi2 = np.random.uniform(-1.0, 1.0)

        xi[i] = [xi1, xi2, xi3]

    return xi


# ============================================================
# 2. Collapsed prismatic vertex basis (BOOK: Fig. D.3(b))
# ============================================================

def prism_vertex_basis_book(xi, eps=1e-12):
    xi1, xi2, xi3 = xi

    denom = 1.0 - xi3
    if abs(denom) < eps:
        denom = eps

    eta_bar = 2.0 * (1.0 + xi1) / denom - 1.0

    L0 = 0.5 * (1.0 - eta_bar)
    L1 = 0.5 * (1.0 + eta_bar)

    M0 = 0.5 * (1.0 - xi2)
    M1 = 0.5 * (1.0 + xi2)

    N0 = 0.5 * (1.0 - xi3)
    N1 = 0.5 * (1.0 + xi3)

    # Vertex order: [A, B, C, D, E, F]
    NA = L0 * M0 * N0
    NB = L1 * M0 * N0
    NC = L0 * M1 * N0

    ND = L1 * M1 * N0   # opposite of A
    NE = M0 * N1
    NF = M1 * N1

    return np.array([NA, NB, NC, ND, NE, NF])


# ============================================================
# 3. Mapping: reference cube -> physical prism
# ============================================================

def map_prism_from_cube(xi, X6):
    N = prism_vertex_basis_book(xi)
    x = (N[:, None] * X6).sum(axis=0)
    return x


# ============================================================
# 4. Geometry helpers
# ============================================================

def prism_edges_book():
    """
    Prism edges following Fig. D.3(b):

      Triangles:
        - A B E
        - D C F

      Connecting edges:
        - A D
        - B C
        - E F
    """
    A, B, C, D, E, F = 0, 1, 2, 3, 4, 5
    return [
        (A, B), (B, E), (E, A),
        (D, C), (C, F), (F, D),
        (A, D), (B, C), (E, F),
    ]


def plot_prism_wire(ax, X6):
    for i, j in prism_edges_book():
        ax.plot([X6[i, 0], X6[j, 0]],
                [X6[i, 1], X6[j, 1]],
                [X6[i, 2], X6[j, 2]],
                color="k", linewidth=1.5)


def label_vertices(ax, X6, labels, offset=(0.03, 0.03, 0.03)):
    for p, name in zip(X6, labels):
        ax.text(p[0] + offset[0],
                p[1] + offset[1],
                p[2] + offset[2],
                name, fontsize=12, fontweight="bold")


# ============================================================
# 5. Main
# ============================================================

def main():

    # --------------------------------------------------------
    # Physical prism (BOOK-consistent geometry)
    # --------------------------------------------------------
    XA = np.array([0.0, 0.0, 0.0])
    XB = np.array([2.0, 0.2, 0.3])
    XE = np.array([0.5, 1.4, 0.1])

    shift = np.array([0.3, 0.4, 1.6])

    XD = XA + shift
    XC = XB + shift
    XF = XE + shift

    # Vertex order: [A, B, C, D, E, F]
    X6 = np.vstack([XA, XB, XC, XD, XE, XF])

    # --------------------------------------------------------
    # Random points in reference prism (xi-coordinates)
    # --------------------------------------------------------
    npts = 600
    xi_pts = random_points_reference_prism(npts)

    # Map to physical
    x_phys = np.array([map_prism_from_cube(xi, X6) for xi in xi_pts])

    # --------------------------------------------------------
    # Plot 1: Reference cube + reference prism points
    # --------------------------------------------------------
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection="3d")

    ax1.scatter(xi_pts[:, 0], xi_pts[:, 1], xi_pts[:, 2],
                s=10, alpha=0.8)

    ax1.set_title("Reference prism embedded in reference cube")
    ax1.set_xlim([-1, 1])
    ax1.set_ylim([-1, 1])
    ax1.set_zlim([-1, 1])
    ax1.set_xlabel(r"$\xi_1$")
    ax1.set_ylabel(r"$\xi_2$")
    ax1.set_zlabel(r"$\xi_3$")

    # --------------------------------------------------------
    # Plot 2: Physical prism + mapped points
    # --------------------------------------------------------
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection="3d")

    plot_prism_wire(ax2, X6)
    ax2.scatter(X6[:, 0], X6[:, 1], X6[:, 2],
                s=70, color="red", label="Vertices")

    label_vertices(ax2, X6, ["A", "B", "C", "D", "E", "F"])

    ax2.scatter(x_phys[:, 0], x_phys[:, 1], x_phys[:, 2],
                s=10, alpha=0.7, label="Mapped points")

    ax2.set_title("Physical prism (collapsed mapping from ref cube)")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_zlabel("z")
    ax2.legend()

    plt.show()


if __name__ == "__main__":
    main()
