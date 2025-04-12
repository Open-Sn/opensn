import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.art3d as art3d
import numpy as np
import math
import sys


def plot_sldfe_quadrature(fileBase):

    print(f"Processing file base: {fileBase}")

    vertPath = fileBase + "_verts.csv"
    cellPath = fileBase + "_cells.csv"
    pointPath = fileBase + "_points.csv"

    # Read vertices
    verts = []
    with open(vertPath, 'r') as verts_file:
        next(verts_file)  # Skip header line
        for line in verts_file:
            # Strip whitespace and split the line by commas, converting each element to a float
            values = [float(v) for v in line.strip().split(',')]
            verts.append(np.array(values))

    # Read cells
    with open(cellPath, 'r') as cells_file:
        next(cells_file)  # Skip header line
        cells = [
            [int(word) for word in line.strip().split(",")[:-1]]
            for line in cells_file
        ]

    # Read points
    weightsum = 0.0
    with open(pointPath, 'r') as points_file:
        next(points_file)  # Skip header line
        points = []
        for line in points_file:
            # Convert each value in the line to a float
            point = [float(val) for val in line.strip().split(",")]
            points.append(point)
            weightsum += point[3]

    print("Weightsum check: ", weightsum, weightsum / 4 / math.pi)

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


if __name__ == "__main__":
    # When run as a script, check command-line arguments or prompt the user
    if len(sys.argv) > 1:
        fileBase = sys.argv[1]
    else:
        fileBase = input("Enter the file base name: ")

    plot_sldfe_quadrature(fileBase)
