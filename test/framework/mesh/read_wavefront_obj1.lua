-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.51187 and 1.42458e-03
num_procs = 4
--Unstructured mesh




--############################################### Check num_procs
if (check_num_procs==nil and number_of_processes ~= num_procs) then
    log.Log(LOG_0ERROR, "Incorrect amount of processors. " ..
                        "Expected "..tostring(num_procs)..
                        ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
meshgen1 = mesh.MeshGenerator.Create
({
  inputs =
  {
    mesh.FromFileMeshGenerator.Create
    ({
      filename = "reactor_pin_mesh.obj"
    })
  }
})
mesh.MeshGenerator.Execute(meshgen1)
--############################################### Exports
if master_export == nil then
    mesh.ExportToVTK("ZObjMesh")
end
