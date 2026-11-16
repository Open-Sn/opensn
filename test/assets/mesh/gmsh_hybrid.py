#!/usr/bin/env python3
"""Generate gmsh_hybrid_v22.msh, gmsh_hybrid_v41.msh, and gmsh_hybrid_v41_binary.msh.

The mesh holds every 3D element type that Gmsh writes in first order: hexahedra in the box
[0, 1]^3, tetrahedra in the box [1, 2] x [0, 1] x [0, 1], with the pyramids that Gmsh inserts
where they meet the quadrilateral faces of the hexahedra, and prisms in [1, 2] x [0, 1] x [1, 1.5],
extruded in two layers from the top of the tetrahedral box. The volume is 2.5 and the boundary
area is 12. All volumes are physical volume 1 and all boundary surfaces are physical surface 2,
named "outer". Run with the gmsh Python module; the supplied files were made with Gmsh 4.15.2.
"""

import gmsh

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.option.setNumber("General.NumThreads", 1)
gmsh.model.add("gmsh_hybrid")
occ = gmsh.model.occ
hex_box = occ.addBox(0, 0, 0, 1, 1, 1)
tet_box = occ.addBox(1, 0, 0, 1, 1, 1)
occ.fragment([(3, hex_box)], [(3, tet_box)])
occ.synchronize()

# Structured, recombined hexahedra in the first box.
for _, surface in gmsh.model.getBoundary([(3, hex_box)], combined=False, oriented=False):
    for _, curve in gmsh.model.getBoundary([(2, surface)], combined=False, oriented=False):
        gmsh.model.mesh.setTransfiniteCurve(curve, 4)
    gmsh.model.mesh.setTransfiniteSurface(surface)
    gmsh.model.mesh.setRecombine(2, surface)
gmsh.model.mesh.setTransfiniteVolume(hex_box)
gmsh.model.mesh.setRecombine(3, hex_box)

# Prisms extruded from the top of the tetrahedral box.
top = [tag for _, tag in gmsh.model.getBoundary([(3, tet_box)], combined=False, oriented=False)
       if abs(gmsh.model.occ.getCenterOfMass(2, tag)[2] - 1.0) < 1e-9]
occ.extrude([(2, top[0])], 0, 0, 0.5, numElements=[2], recombine=True)
occ.synchronize()

volumes = [tag for _, tag in gmsh.model.getEntities(3)]
gmsh.model.addPhysicalGroup(3, volumes, 1)
boundary = gmsh.model.getBoundary([(3, v) for v in volumes], combined=True, oriented=False)
gmsh.model.addPhysicalGroup(2, sorted(abs(tag) for _, tag in boundary), 2, "outer")
gmsh.option.setNumber("Mesh.MeshSizeMax", 0.35)
gmsh.model.mesh.generate(3)

types, tags, _ = gmsh.model.mesh.getElements(3)
print({int(t): len(x) for t, x in zip(types, tags)})
for name, version, binary in (("gmsh_hybrid_v22.msh", 2.2, 0), ("gmsh_hybrid_v41.msh", 4.1, 0),
                              ("gmsh_hybrid_v41_binary.msh", 4.1, 1)):
    gmsh.option.setNumber("Mesh.MshFileVersion", version)
    gmsh.option.setNumber("Mesh.Binary", binary)
    gmsh.write(name)
gmsh.finalize()
