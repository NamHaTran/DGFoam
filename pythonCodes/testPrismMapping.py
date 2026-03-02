# ---------------------------------------------------------------
# Prism visualization with equal axis scaling AND equal tick spacing
# - Single plot
# - No manual color specification
# ---------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Prism vertices
# ---------------------------------------------------------------
owner_vertices = np.array([
    [17.0792, 74.1771, 24.6361],
    [18.4082, 75.434 , 23.8272],
    [17.5629, 72.7914, 23.2775],
    [15.7531, 74.8402, 23.4876],
    [17.082 , 76.097 , 22.6787],
    [16.2368, 73.4544, 22.129 ]
])

neighbor_vertices = np.array([
    [17.5629, 72.7914, 23.2775],
    [18.4082, 75.434 , 23.8272],
    [19.4915, 74.8623, 22.2462],
    [16.2368, 73.4544, 22.129 ],
    [17.082 , 76.097 , 22.6787],
    [18.1654, 75.5254, 21.0977]
])

# Prism edges
edges = [
    (0,1),(1,2),(2,0),
    (3,4),(4,5),(5,3),
    (0,3),(1,4),(2,5)
]

# Gauss points
owner_phys = np.array([
    [17.9493, 75.0156, 23.4683],
    [17.4613, 73.4899, 23.151 ],
    [17.1836, 75.3985, 22.8053],
    [16.6956, 73.8728, 22.4879]
])

neighbor_phys = np.array([
    [17.9493, 75.0156, 23.4683],
    [17.1836, 75.3985, 22.8053],
    [17.4613, 73.4899, 23.151 ],
    [16.6956, 73.8728, 22.4879]
])

# ---------------------------------------------------------------
# Plot
# ---------------------------------------------------------------
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Draw edges
for verts in [owner_vertices, neighbor_vertices]:
    for i, j in edges:
        ax.plot(
            [verts[i,0], verts[j,0]],
            [verts[i,1], verts[j,1]],
            [verts[i,2], verts[j,2]]
        )

# Draw Gauss points
ax.scatter(owner_phys[:,0], owner_phys[:,1], owner_phys[:,2])
ax.scatter(neighbor_phys[:,0], neighbor_phys[:,1], neighbor_phys[:,2])

# ---------------------------------------------------------------
# Equal axis scaling
# ---------------------------------------------------------------
all_pts = np.vstack((owner_vertices, neighbor_vertices))
mins = all_pts.min(axis=0)
maxs = all_pts.max(axis=0)

ranges = maxs - mins
max_range = ranges.max()

centers = (maxs + mins) / 2.0

ax.set_xlim(centers[0] - max_range/2, centers[0] + max_range/2)
ax.set_ylim(centers[1] - max_range/2, centers[1] + max_range/2)
ax.set_zlim(centers[2] - max_range/2, centers[2] + max_range/2)

# ---------------------------------------------------------------
# Equal tick spacing
# ---------------------------------------------------------------
n_ticks = 6
tick_vals = np.linspace(-max_range/2, max_range/2, n_ticks)

ax.set_xticks(centers[0] + tick_vals)
ax.set_yticks(centers[1] + tick_vals)
ax.set_zticks(centers[2] + tick_vals)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

plt.show()
