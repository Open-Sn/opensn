#include "opensn/framework/mesh/MeshModifiers/SnapToPlaneMeshModifier.h"

#include "opensn/framework/ChiObjectFactory.h"

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/mpi/chi_mpi.h"
#include "opensn/framework/logging/chi_log.h"

#include "opensn/framework/mesh/MeshHandler/chi_meshhandler.h"
#include "opensn/framework/mesh/MeshContinuum/chi_meshcontinuum.h"

namespace chi_mesh
{

RegisterChiObject(chi_mesh, SnapToPlaneMeshModifier);

chi::InputParameters
SnapToPlaneMeshModifier::GetInputParameters()
{
  chi::InputParameters params = ChiObject::GetInputParameters();

  params.SetGeneralDescription(
    "Modifier that will snap nodes, that are within a tolerated distance of a "
    "plane, to the plane. This modifier is useful for straitening a boundary "
    "edge or aligning vertices that misaligned during meshing.");
  params.SetDocGroup("DocMeshModifiers");

  params.AddRequiredParameterArray("normal",
                                   "The normal of the plane to which the nodes are to be snapped.");

  params.AddRequiredParameterArray("point", "The anchor point of the plane.");

  params.AddOptionalParameter(
    "boundaries_only", true, "If set to true, only boundary nodes will be snapped.");

  params.AddOptionalParameter(
    "check_face_alignment",
    false,
    "If set, only faces that match the plane normal (within tolerance) will "
    "have their nodes snapped.");

  params.AddOptionalParameter(
    "tolerance", 1.0e-5, "Tolerance per dimension within which a face/edge will be aligned");

  return params;
}

SnapToPlaneMeshModifier::SnapToPlaneMeshModifier(const chi::InputParameters& params)
  : MeshModifier(params),
    normal_(params.GetParamVectorValue<double>("normal")),
    point_(params.GetParamVectorValue<double>("point")),
    boundary_nodes_only_(params.GetParamValue<bool>("boundaries_only")),
    check_face_alignment_(params.GetParamValue<bool>("check_face_alignment")),
    tol_(params.GetParamValue<double>("tolerance"))
{
  ChiLogicalErrorIf(Chi::mpi.process_count != 1, "Cannot only be used in serial");
}

void
SnapToPlaneMeshModifier::Apply()
{
  auto& grid = *chi_mesh::GetCurrentHandler().GetGrid();

  std::vector<uint64_t> snapped_vertex_ids;
  std::set<uint64_t> cell_ids_modified;

  if (check_face_alignment_)
    for (const auto& cell : grid.local_cells)
    {
      for (const auto& face : cell.faces_)
      {
        if (boundary_nodes_only_ and face.has_neighbor_) continue;

        bool matches = true;
        for (size_t d = 0; d < 3; ++d)
          if (std::fabs(face.normal_[d] - normal_[d]) >= tol_) matches = false;

        if (matches)
          for (uint64_t vid : face.vertex_ids_)
          {
            const double d = (grid.vertices[vid] - point_).Dot(normal_);

            if (std::fabs(d) < tol_)
            {
              grid.vertices[vid] -= d * normal_;
              snapped_vertex_ids.push_back(vid);
              cell_ids_modified.insert(cell.local_id_);
            }
          }
      } // for face
    }   // for cell
  else
    for (const auto& cell : grid.local_cells)
    {
      for (const uint64_t vid : cell.vertex_ids_)
      {
        const double d = (grid.vertices[vid] - point_).Dot(normal_);

        if (std::fabs(d) < tol_)
        {
          grid.vertices[vid] -= d * normal_;
          snapped_vertex_ids.push_back(vid);
          cell_ids_modified.insert(cell.local_id_);
        }
      }
    } // for cell

  // Modifying cells
  for (const uint64_t cell_local_id : cell_ids_modified)
  {
    grid.local_cells[cell_local_id].RecomputeCentroidsAndNormals(grid);
    // for (const auto& face : grid.local_cells[cell_local_id_].faces_)
    //   chi::log.Log() << face.normal_.PrintStr();
  }

  Chi::log.Log0Verbose1() << "Number of cells modified " << cell_ids_modified.size();
}

} // namespace chi_mesh
