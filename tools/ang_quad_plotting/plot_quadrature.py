import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
from pathlib import Path
import numpy as np
import math

# Declare File Base to read
fileBase = "TestQuad1_"

currentDirect = Path(__file__).parent
vertPath = currentDirect / (fileBase + "verts.csv")
cellPath = currentDirect / (fileBase + "cells.csv")
pointPath = currentDirect / (fileBase + "points.csv")

# Read vertices
verts = []
verts_file = open(vertPath)
verts_file.readline()
for line in verts_file:
    words = line.split(",")
    verts.append(np.array([float(words[0]), float(words[1]), float(words[2])]))
verts_file.close()

# Read cells
cells = []
cells_file = open(cellPath)
cells_file.readline()
for line in cells_file:
    words = line.split(",")
    words.pop(-1)
    cell = []
    for word in words:
        cell.append(int(word))
    cells.append(cell)
cells_file.close()

# Read points
points = []
weightsum = 0.0
points_file = open(pointPath)
points_file.readline()
for line in points_file:
    words = line.split(",")
    point = []
    for word in words:
        point.append(float(word))
    points.append(point)
    weightsum += point[3]
points_file.close()

print("Weightsum check: ", weightsum, weightsum / 4 / math.pi)

points_array = np.array(points)

# Generate polygons
patches = []

for cell in cells:

    vertex_list = []

    for index in cell:
        vertex_list.append(verts[index])

    polygon = art3d.Poly3DCollection([vertex_list])
    polygon.set_color([1.0, 1.0, 1.0, 1.0])
    polygon.set_edgecolor([0.0, 0.0, 0.0, 1.0])
    patches.append(polygon)

# Plot polygons
fig = plt.figure(figsize=(10, 8.5))
ax = fig.add_subplot(111, projection='3d')

ax.view_init(20, 45)
limit = 1

for poly in patches:
    ax.add_collection3d(poly)

if limit == 8:
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.0])
    ax.set_zlim([0.0, 1.0])
else:
    ax.set_xlim([-1.0, 1.0])
    ax.set_ylim([-1.0, 1.0])
    ax.set_zlim([-1.0, 1.0])

ax.margins(0.5)
ax.set_xlabel(r"$\mu$")
ax.set_ylabel(r"$\eta$")
ax.set_zlabel(r"$\xi$")
plt.show()
