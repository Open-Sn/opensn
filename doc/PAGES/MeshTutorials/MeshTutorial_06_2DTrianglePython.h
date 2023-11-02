/** \page MeshTutorial_06 Mesh Tutorial 6: A 2D Unstructured Mesh

 Here we provide an example of a 2D unstructured mesh generated using the
 2D ```Triangle``` mesh generator, with is available in Python via the ```MeshPy``` module.
 We export the resulting mesh using ```meshio``` in both ```.obj``` and ```.vtu```
 formats. Note that only the export as VTU contains material IDs.

\code
import numpy as np
import meshpy.triangle as triangle
import meshio

n_circumferential_points =20
radius = 1.0
width = 3.0

points = []
# Define the vertices of the square domain (centered on 0)
points = [(-width/2,-width/2), (width/2,-width/2),  (width/2,width/2),  (-width/2,width/2),]
# Create a list of points for the circle (centered on 0)
for i in range(n_circumferential_points):
    angle = i * 2 * np.pi / n_circumferential_points
    points.append((radius * np.cos(angle), radius * np.sin(angle)))

# Define the segments of the square domain
segments = [(0,1), (1,2), (2,3), (3,0)]
markers = [1,2,3,4]
Nbeg = len(markers)
# Define the segments for the circle
for i in range(Nbeg,n_circumferential_points+Nbeg):
    if i+1 == n_circumferential_points+Nbeg:
        ip1 = Nbeg
    else:
        ip1 = i+1
    segments.extend([(i,ip1)])
    markers.extend([0])

# Create a mesh info object
mesh_info = triangle.MeshInfo()
# Add points
mesh_info.set_points(points)
# Add facets
mesh_info.set_facets(segments, facet_markers=markers)
# Add regions
mesh_info.regions.resize(2)
# inside circle region
mesh_info.regions[0] = [0.0, 0.0, 1, 0.1] # 2d pt + ID + max_vol
# outside circle region
eps = 1e-2
mesh_info.regions[1] = [-width/2+eps,-width/2+eps, 2, 0.1] # 2d pt + ID + max_vol

# Refine the mesh
mesh = triangle.build(mesh_info, verbose=True, max_volume=0.02, min_angle=25,
                      attributes=True, generate_faces=True)

# Add a third dimension to the points
points_3d = [(x, y, 0) for x, y in mesh.points]

directory = ""
filename = "pin_01"
meshio.write_points_cells(directory + filename + ".obj",
                          points_3d,
                          {"triangle": mesh.elements})

# Save the mesh to an .vtk file
meshio.write_points_cells(directory + filename + ".vtu",
                          points_3d,
                          {"triangle": mesh.elements},
                          cell_data={"attribute": [mesh.element_attributes]})

 \endcode


*/