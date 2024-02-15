/**
\addtogroup chi_mesh__RCCLogicalVolume

## Additional Example
\code
MeshHandlerCreate()

mesh={}
N=40
L=5
xmin = -L/2
dx = L/N
for i=1,(N+1) do
k=i-1
mesh[i] = xmin + k*dx
end
MeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,mesh)
VolumeMesherExecute();

lv1 = chi_mesh.RCCLogicalVolume.Create({r = 1.3, x0=L/2, y0=L/2, z0 = -1.0, vz
= 2.0}) VolumeMesherSetProperty(MATID_FROMLOGICAL, lv1, 1)

lv2 = chi_mesh.RCCLogicalVolume.Create({r = 1.3,
x0=-0.8, y0=-0.8, z0=-1.5,
vx=1.0, vy=1.0, vz=3.0})
VolumeMesherSetProperty(MATID_FROMLOGICAL, lv2, 2)

mesh.ExportToVTK("lv_rcc_test1")
\endcode

\image html framework/chi_mesh/LogicalVolume/lv_rcc_test1.png width=500px
*/
