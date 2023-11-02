#include "lbs_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

#include <fstream>
#include <cstring>

//###################################################################
/**Writes the groupset's angular fluxes to file.*/
void
lbs::LBSSolver::WriteGroupsetAngularFluxes(const LBSGroupset& groupset,
                                           const std::string& file_base)
{
  std::string file_name = file_base + std::to_string(Chi::mpi.location_id) + ".data";

  //============================================= Open file
  std::ofstream file(file_name,
                     std::ofstream::binary |  // binary file
                       std::ofstream::out |   // no accidental reading
                       std::ofstream::trunc); // clear file contents when opened

  //============================================= Check file is open
  if (not file.is_open())
  {
    Chi::log.LogAllWarning() << __FUNCTION__ << "Failed to open " << file_name;
    return;
  }

  //============================================= Write header
  std::string header_info = "Chi-Tech LinearBoltzmann::Groupset angular flux file\n"
                            "Header size: 320 bytes\n"
                            "Structure(type-info):\n"
                            "size_t-num_local_nodes\n"
                            "size_t-num_angles\n"
                            "size_t-num_groups\n"
                            "size_t-num_records\n"
                            "Each record:\n"
                            "size_t-cell_global_id\n"
                            "unsigned int-node_number\n"
                            "unsigned int-angle_num_\n"
                            "unsigned int-group_num\n"
                            "double-angular_flux\n";

  int header_size = (int)header_info.length();

  char header_bytes[320];
  memset(header_bytes, '-', 320);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, 319));
  header_bytes[319] = '\0';

  file << header_bytes;

  //============================================= Get relevant items
  auto NODES_ONLY = chi_math::UnknownManager::GetUnitaryUnknownManager();

  size_t num_local_nodes = discretization_->GetNumLocalDOFs(NODES_ONLY);
  size_t num_angles = groupset.quadrature_->abscissae_.size();
  size_t num_groups = groupset.groups_.size();
  size_t num_local_dofs = psi_new_local_[groupset.id_].size();
  auto dof_handler = groupset.psi_uk_man_;

  //============================================= Write num_ quantities
  file.write((char*)&num_local_nodes, sizeof(size_t));
  file.write((char*)&num_angles, sizeof(size_t));
  file.write((char*)&num_groups, sizeof(size_t));
  file.write((char*)&num_local_dofs, sizeof(size_t));

  auto& sdm = discretization_;

  //============================================= Write per dof data
  size_t dof_count = 0;
  for (const auto& cell : grid_ptr_->local_cells)
  {
    const size_t num_nodes = sdm->GetCellNumNodes(cell);
    for (unsigned int i = 0; i < num_nodes; ++i)
      for (unsigned int n = 0; n < num_angles; ++n)
        for (unsigned int g = 0; g < num_groups; ++g)
        {
          if (++dof_count > num_local_dofs) goto close_file;

          uint64_t dof_map = sdm->MapDOFLocal(cell, i, dof_handler, n, g);
          double value = psi_new_local_[groupset.id_][dof_map];

          file.write((char*)&cell.global_id_, sizeof(size_t));
          file.write((char*)&i, sizeof(unsigned int));
          file.write((char*)&n, sizeof(unsigned int));
          file.write((char*)&g, sizeof(unsigned int));
          file.write((char*)&value, sizeof(double));
        }
  }

//============================================= Clean-up
close_file:
  file.close();
}

//###################################################################
/**Prints the groupset's angular fluxes to file.*/
void
lbs::LBSSolver::ReadGroupsetAngularFluxes(LBSGroupset& groupset, const std::string& file_base)
{
  std::string file_name = file_base + std::to_string(Chi::mpi.location_id) + ".data";

  //============================================= Open file
  std::ifstream file(file_name,
                     std::ofstream::binary | // binary file
                       std::ofstream::in);   // no accidental writing

  //============================================= Check file is open
  if (not file.is_open())
  {
    Chi::log.LogAllWarning() << __FUNCTION__ << "Failed to open " << file_name;
    return;
  }

  //============================================= Get relevant items
  auto NODES_ONLY = chi_math::UnknownManager::GetUnitaryUnknownManager();

  size_t num_local_nodes = discretization_->GetNumLocalDOFs(NODES_ONLY);
  size_t num_angles = groupset.quadrature_->abscissae_.size();
  size_t num_groups = groupset.groups_.size();
  size_t num_local_dofs = psi_new_local_[groupset.id_].size();
  std::vector<double>& psi = psi_new_local_[groupset.id_];
  auto dof_handler = groupset.psi_uk_man_;

  size_t file_num_local_nodes;
  size_t file_num_angles;
  size_t file_num_groups;
  size_t file_num_local_dofs;

  //============================================= Read header
  Chi::log.Log() << "Reading angular flux file " << file_name;
  char header_bytes[320];
  header_bytes[319] = '\0';
  file.read(header_bytes, 319);

  file.read((char*)&file_num_local_nodes, sizeof(size_t));
  file.read((char*)&file_num_angles, sizeof(size_t));
  file.read((char*)&file_num_groups, sizeof(size_t));
  file.read((char*)&file_num_local_dofs, sizeof(size_t));

  //============================================= Check compatibility
  if (file_num_local_nodes != num_local_nodes or file_num_angles != num_angles or
      file_num_groups != num_groups or file_num_local_dofs != num_local_dofs)
  {
    std::stringstream outstr;
    outstr << "num_local_nodes: " << file_num_local_nodes << "\n";
    outstr << "num_angles     : " << file_num_angles << "\n";
    outstr << "num_groups     : " << file_num_groups << "\n";
    outstr << "num_local_dofs : " << file_num_local_dofs << "\n";
    Chi::log.LogAll() << "Incompatible DOF data found in file " << file_name << "\n"
                      << outstr.str();
    file.close();
    return;
  }

  auto& sdm = discretization_;

  //============================================= Commit to reading the file
  psi.reserve(file_num_local_dofs);
  std::set<uint64_t> cells_touched;
  for (size_t dof = 0; dof < file_num_local_dofs; ++dof)
  {
    uint64_t cell_global_id;
    unsigned int node;
    unsigned int angle_num;
    unsigned int group;
    double psi_value;

    file.read((char*)&cell_global_id, sizeof(uint64_t));
    file.read((char*)&node, sizeof(unsigned int));
    file.read((char*)&angle_num, sizeof(unsigned int));
    file.read((char*)&group, sizeof(unsigned int));
    file.read((char*)&psi_value, sizeof(double));

    cells_touched.insert(cell_global_id);

    const auto& cell = grid_ptr_->cells[cell_global_id];

    size_t imap = sdm->MapDOFLocal(cell, node, dof_handler, angle_num, group);

    psi[imap] = psi_value;
  }

  Chi::log.LogAll() << "Number of cells read: " << cells_touched.size();

  //============================================= Clean-up
  file.close();
}