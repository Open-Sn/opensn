import os

# Define the vertices, faces, and normals for a cube again
vertices = [
    [-0.5, -0.5, -0.5], [0.5, -0.5, -0.5],
    [0.5, 0.5, -0.5], [-0.5, 0.5, -0.5],
    [-0.5, -0.5, 0.5], [0.5, -0.5, 0.5],
    [0.5, 0.5, 0.5], [-0.5, 0.5, 0.5]
]

faces = [
    [1, 2, 3, 4], [5, 8, 7, 6],
    [1, 5, 6, 2], [2, 6, 7, 3],
    [3, 7, 8, 4], [5, 1, 4, 8]
]

normals = [
    [0, 0, -1], [0, 0, 1],
    [0, -1, 0], [0, 1, 0],
    [-1, 0, 0], [1, 0, 0]
]

# Path for the new .obj file including normals
obj_file_path_with_normals = "./cube_with_normals.obj"

# Create or overwrite the .obj file to include normals
with open(obj_file_path_with_normals, "w") as file:
    file.write("# Cube with Normals\n")
    for v in vertices:
        file.write(f"v {' '.join(map(str, v))}\n")
    for n in normals:
        file.write(f"vn {' '.join(map(str, n))}\n")
    # Each face needs to be updated to use the correct normal index
    # Assuming each side of the cube has a unique normal, we map face index to normal index
    for i, f in enumerate(faces):
        normal_index = i // 2 + 1  # Simplified mapping for cube faces to normals
        file.write(f"f {' '.join([f'{v}//{normal_index}' for v in f])}\n")

# Check if the file has been created and include normals
os.path.exists(obj_file_path_with_normals), obj_file_path_with_normals
