#include "lbsadj_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

#include "math/chi_math_vectorNX.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "lbs_adjoint.h"

#include <fstream>

//###################################################################
/**Exports an importance map in binary format.*/
void
lbs::DiscreteOrdinatesAdjointSolver::ExportImportanceMap(const std::string& file_name)
{
  const std::string fname = __FUNCTION__;

  //============================================= Determine cell averaged
  //                                              importance map
  std::set<int> set_group_numbers;
  for (const auto& groupset : groupsets_)
    for (const auto& group : groupset.groups_)
      set_group_numbers.insert(group.id_);

  const auto& m_to_ell_em_map = groupsets_.front().quadrature_->GetMomentToHarmonicsIndexMap();

  typedef chi_math::VectorN<4> Arr4; // phi, J_x, J_y, J_z
  typedef std::vector<Arr4> MGVec4;
  typedef std::vector<MGVec4> VecOfMGVec4;
  const size_t num_groups = set_group_numbers.size();
  const size_t num_cells = grid_ptr_->local_cells.size();

  VecOfMGVec4 cell_avg_p1_moments(num_cells, MGVec4(num_groups));
  {

    for (const auto& cell : grid_ptr_->local_cells)
    {
      const auto& cell_view = cell_transport_views_[cell.local_id_];
      const int num_nodes = cell_view.NumNodes();
      const auto& fe_values = unit_cell_matrices_[cell.local_id_];

      VecOfMGVec4 nodal_p1_moments(num_nodes);
      for (int i = 0; i < num_nodes; ++i)
      {
        //==================================== Get multigroup p1_moments
        MGVec4 p1_moments(set_group_numbers.size(), VecDbl{0.0, 0.0, 0.0, 0.0});
        for (int m = 0; m < std::max(static_cast<int>(num_moments_), 4); ++m)
        {
          const auto& ell = m_to_ell_em_map[m].ell;
          const auto& em = m_to_ell_em_map[m].m;

          size_t dof_map_g0 = cell_view.MapDOF(i, m, 0); // unknown map

          for (int g : set_group_numbers)
          {
            if (ell == 0 and em == 0) p1_moments[g](0) = std::fabs(phi_old_local_[dof_map_g0 + g]);
            if (ell == 1 and em == 1) p1_moments[g](1) = phi_old_local_[dof_map_g0 + g];
            if (ell == 1 and em == -1) p1_moments[g](2) = phi_old_local_[dof_map_g0 + g];
            if (ell == 1 and em == 0) p1_moments[g](3) = phi_old_local_[dof_map_g0 + g];
          } // for g
        }   // for m

        nodal_p1_moments[i] = std::move(p1_moments);
      } // for node i

      //=========================================== Determine nodal average
      //                                            p1_moments
      for (int g : set_group_numbers)
      {
        chi_math::VectorN<4> cell_p1_avg(VecDbl{0.0, 0.0, 0.0, 0.0});

        double volume_total = 0.0;
        for (int i = 0; i < num_nodes; ++i)
        {
          double IntV_shapeI = fe_values.Vi_vectors[i];
          cell_p1_avg += nodal_p1_moments[i][g] * IntV_shapeI;
          volume_total += IntV_shapeI;
        } // for node i
        cell_p1_avg /= volume_total;

        cell_avg_p1_moments[cell.local_id_][g] = cell_p1_avg;
      } // for g
    }   // for cell
  }

  //============================================= Determine cell-based
  //                                              exponential-representations
  typedef std::pair<double, double> ABCoeffsPair;
  typedef std::vector<ABCoeffsPair> VecOfABCoeffsPair;
  typedef std::vector<VecOfABCoeffsPair> ExpReps;
  ExpReps cell_exp_reps(num_cells, VecOfABCoeffsPair(num_groups, {0.0, 0.0}));
  {
    for (const auto& cell : grid_ptr_->local_cells)
    {
      for (int g : set_group_numbers)
      {
        const auto& p1_moments = cell_avg_p1_moments[cell.local_id_][g];

        auto a_b = MakeExpRepFromP1({p1_moments[0], p1_moments[1], p1_moments[2], p1_moments[3]});

        cell_exp_reps[cell.local_id_][g] = std::make_pair(a_b[0], a_b[1]);
      } // for g
    }   // for cell
  }

  Chi::log.Log() << "Exporting importance map to binary file " << file_name;

  const auto locJ_io_flags = std::ofstream::binary | std::ofstream::out;
  const auto loc0_io_flags = locJ_io_flags | std::ofstream::trunc;
  const bool is_home = (Chi::mpi.location_id == 0);

  //======================================== Build header
  std::string header_info = "Chi-Tech LinearBoltzmann: Importance map file\n"
                            "Header size: 500 bytes\n"
                            "Structure(type-info):\n"
                            "uint64_t num_global_cells\n"
                            "uint64_t num_groups\n"
                            "uint64_t num_records\n"
                            "Each record:\n"
                            "  uint64_t     cell_global_id\n"
                            "  unsigned int group_num\n"
                            "  double       phi\n"
                            "  double       J_x\n"
                            "  double       J_y\n"
                            "  double       J_z\n"
                            "  double       a_coefficient\n"
                            "  double       b_coefficient\n";

  int header_size = (int)header_info.length();

  char header_bytes[400];
  memset(header_bytes, '-', 400);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, 399));
  header_bytes[399] = '\0';

  //================================================== Process each location
  uint64_t num_global_cells = grid_ptr_->GetGlobalNumberOfCells();
  for (int locationJ = 0; locationJ < Chi::mpi.process_count; ++locationJ)
  {
    Chi::log.LogAll() << "  Barrier at " << locationJ;
    Chi::mpi.Barrier();
    if (Chi::mpi.location_id != locationJ) continue;

    Chi::log.LogAll() << "  Location " << locationJ << " appending data.";

    std::ofstream file(file_name, is_home ? loc0_io_flags : locJ_io_flags);

    if (not file.is_open())
    {
      std::stringstream outstr;

      outstr << fname << ": Location " << Chi::mpi.location_id << ", failed to open file "
             << file_name;
      throw std::logic_error(outstr.str());
    }

    if (is_home)
    {
      file << header_bytes;
      uint64_t num_groups_t = groups_.size();
      uint64_t num_records = num_global_cells * num_groups;

      file.write((char*)&num_global_cells, sizeof(uint64_t));
      file.write((char*)&num_groups_t, sizeof(uint64_t));
      file.write((char*)&num_records, sizeof(uint64_t));
    }

    for (const auto& cell : grid_ptr_->local_cells)
    {
      MGVec4 p1_moments = cell_avg_p1_moments[cell.local_id_];
      VecOfABCoeffsPair exp_rep = cell_exp_reps[cell.local_id_];

      auto cell_global_id = static_cast<uint64_t>(cell.global_id_);
      for (int group : set_group_numbers)
      {
        auto g = static_cast<unsigned int>(group);
        file.write((char*)&cell_global_id, sizeof(uint64_t));
        file.write((char*)&g, sizeof(unsigned int));

        for (int m = 0; m < 4; ++m)
          file.write((char*)&p1_moments[group](m), sizeof(double));

        file.write((char*)&exp_rep[group].first, sizeof(double));
        file.write((char*)&exp_rep[group].second, sizeof(double));
      } // for g
    }   // for cell

    file.close();
  } // for location

  Chi::log.LogAll() << "Done exporting importance map to binary file " << file_name;
  Chi::mpi.Barrier();
}