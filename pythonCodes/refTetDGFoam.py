import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# --------------------------------------------------
# Create figure
# --------------------------------------------------
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, projection='3d')

# --------------------------------------------------
# Reference hexahedron vertices (-1,1)^3
# Index -> Letter mapping: 0-A, 1-B, ..., 7-H
# --------------------------------------------------
letters = ['A', 'B', ' ', 'C', 'E', ' ', ' ', ' ']

vertices = {
    0: (-1, -1, -1),  # A
    1: ( 1, -1, -1),  # B
    #2: ( 1,  1, -1),  # D
    3: (-1,  1, -1),  # C
    4: (-1, -1,  1),  # E
    5: (-1, -1,  1),  # F
    6: (-1, -1,  1),  # G
    7: (-1, -1,  1),  # F
}

edges = [
    (0,1), (3,0),
    (4,5), (5,6), (6,7), (7,4),
    (0,4), (1,5), (3,7), (1,3)
]

# --------------------------------------------------
# Draw edges
# --------------------------------------------------
for i, j in edges:
    xi, yi, zi = vertices[i]
    xj, yj, zj = vertices[j]
    ax.plot([xi, xj], [yi, yj], [zi, zj], color='black', linewidth=1.4)

# --------------------------------------------------
# Draw vertices and labels
# --------------------------------------------------
for idx, (x, y, z) in vertices.items():
    letter = letters[idx]
    ax.scatter(x, y, z, color='blue', s=45)

    label = (
        f"{letter}\n"
        f"({x:+.0f}, {y:+.0f}, {z:+.0f})"
    )

    ax.text(
        x, y, z,
        label,
        color='blue',
        fontsize=9,
        ha='left',
        va='bottom'
    )

# --------------------------------------------------
# Draw reference axes
# --------------------------------------------------
ax.quiver(0, 0, 0, 1.5, 0, 0, color='red', arrow_length_ratio=0.08)
ax.quiver(0, 0, 0, 0, 1.5, 0, color='red', arrow_length_ratio=0.08)
ax.quiver(0, 0, 0, 0, 0, 1.5, color='red', arrow_length_ratio=0.08)

ax.text(1.6, 0, 0, r'$\eta_1$', color='red', fontsize=12)
ax.text(0, 1.6, 0, r'$\eta_2$', color='red', fontsize=12)
ax.text(0, 0, 1.6, r'$\eta_3$', color='red', fontsize=12)

# --------------------------------------------------
# Plot settings
# --------------------------------------------------
ax.set_box_aspect([1, 1, 1])
ax.set_xlim([-2, 2])
ax.set_ylim([-2, 2])
ax.set_zlim([-2, 2])

ax.set_axis_off()
ax.view_init(elev=22, azim=-55)

# --------------------------------------------------
# Save as SVG for Doxygen
# --------------------------------------------------
plt.savefig(
    "/home/nam-ha-tran/OpenFOAM/DGFoam/pythonCodes/images/refTetBasisDomain.svg",
    format="svg",
    bbox_inches="tight"
)
plt.show()
