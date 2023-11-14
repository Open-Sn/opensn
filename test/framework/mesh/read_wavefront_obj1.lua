-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.51187 and 1.42458e-03
num_procs = 4
--Unstructured mesh




--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
meshgen1 = chi_mesh.MeshGenerator.Create
({
  inputs =
  {
    chi_mesh.FromFileMeshGenerator.Create
    ({
      filename = "reactor_pin_mesh.obj"
    })
  }
})
chi_mesh.MeshGenerator.Execute(meshgen1)
--############################################### Exports
if master_export == nil then
    chiMeshHandlerExportMeshToVTK("ZObjMesh")
end
