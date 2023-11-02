#include "sweep_boundaries.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"
#include "math/Quadratures/angular_quadrature_base.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/**Returns a pointer to a heterogeneous flux storage location.*/
double*
chi_mesh::sweep_management::BoundaryIncidentHeterogeneous::HeterogeneousPsiIncoming(
  uint64_t cell_local_id,
  unsigned int face_num,
  unsigned int fi,
  unsigned int angle_num,
  int group_num,
  size_t gs_ss_begin)
{
  if (local_cell_data_.empty())
  {
    Chi::log.LogAllError() << "HeterogeneousPsiIncoming call made to a heterogeneous boundary "
                              "with that information not yet set up.";
    exit(EXIT_FAILURE);
  }

  const size_t dof_offset = num_groups_ * angle_num + group_num;

  return &local_cell_data_[cell_local_id][face_num][fi][dof_offset];
}

//###################################################################
/**Performs the setup for a particular quadrature.*/
void
chi_mesh::sweep_management::BoundaryIncidentHeterogeneous::Setup(
  const chi_mesh::MeshContinuum& grid, const chi_math::AngularQuadrature& quadrature)
{
  const size_t num_local_cells = grid.local_cells.size();
  local_cell_data_.clear();
  local_cell_data_.reserve(num_local_cells);

  std::vector<bool> cell_bndry_flags(num_local_cells, false);
  for (const auto& cell : grid.local_cells)
    for (const auto& face : cell.faces_)
      if (not face.has_neighbor_)
      {
        cell_bndry_flags[cell.local_id_] = true;
        break;
      }

  size_t num_angles = quadrature.omegas_.size();

  typedef std::pair<double, double> PhiTheta;

  std::vector<int> angle_indices;
  std::vector<chi_mesh::Vector3> angle_vectors;
  std::vector<PhiTheta> phi_theta_angles;
  std::vector<int> group_indices;

  angle_indices.reserve(num_angles);
  angle_vectors.reserve(num_angles);
  phi_theta_angles.reserve(num_angles);
  group_indices.reserve(num_groups_);

  int num_angles_int = static_cast<int>(num_angles);
  for (int n = 0; n < num_angles_int; ++n)
    angle_indices.emplace_back(n);
  for (int n = 0; n < num_angles_int; ++n)
    angle_vectors.emplace_back(quadrature.omegas_[n]);
  for (int n = 0; n < num_angles_int; ++n)
  {
    auto& abscissae = quadrature.abscissae_[n];
    double phi = abscissae.phi;
    double theta = abscissae.theta;
    phi_theta_angles.emplace_back(std::make_pair(phi, theta));
  }
  for (int g = 0; g < static_cast<int>(num_groups_); ++g)
    group_indices.emplace_back(g);

  const double eval_time = GetEvaluationTime();

  for (const auto& cell : grid.local_cells)
  {
    if (cell_bndry_flags[cell.local_id_])
    {
      CellData cell_data(cell.faces_.size());

      for (size_t f = 0; f < cell.faces_.size(); ++f)
      {
        auto& face = cell.faces_[f];
        size_t face_num_nodes = face.vertex_ids_.size();
        FaceData face_data;

        if (not face.has_neighbor_ and face.neighbor_id_ == ref_boundary_id_)
        {
          face_data.reserve(face_num_nodes);
          for (size_t i = 0; i < face_num_nodes; ++i)
          {
            std::vector<double> face_node_data =
              boundary_function_->Evaluate(cell.global_id_,
                                           cell.material_id_,
                                           f,
                                           i,
                                           grid.vertices[face.vertex_ids_[i]],
                                           face.normal_,
                                           angle_indices,
                                           angle_vectors,
                                           phi_theta_angles,
                                           group_indices,
                                           eval_time);

            face_data.push_back(std::move(face_node_data));
          } // for face node-i
        }   // bndry face

        cell_data[f] = std::move(face_data);
      } // for face f

      local_cell_data_.push_back(std::move(cell_data));
    } // if bndry cell
    else
      local_cell_data_.emplace_back();

  } // for cell
}