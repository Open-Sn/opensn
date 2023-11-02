--############################################### Setup mesh

meshgen1 = chi_mesh.MeshGenerator.Create
({
  inputs =
  {
    chi_mesh.FromFileMeshGenerator.Create
    ({
      filename = "mesh/2D_c5g7_coarse.msh"
    })
  }
})
chi_mesh.MeshGenerator.Execute(meshgen1)