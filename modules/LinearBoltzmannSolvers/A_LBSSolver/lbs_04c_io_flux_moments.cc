#include "lbs_solver.h"

#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include <fstream>
#include <cstring>
#include <cassert>

//###################################################################
/**Makes a source-moments vector from scattering and fission based
 * on the latest phi-solution.*/
std::vector<double>
lbs::LBSSolver::MakeSourceMomentsFromPhi()
{
  size_t num_local_dofs = discretization_->GetNumLocalDOFs(flux_moments_uk_man_);

  std::vector<double> source_moments(num_local_dofs, 0.0);
  for (auto& groupset : groupsets_)
  {
    active_set_source_function_(groupset,
                                source_moments,
                                PhiOldLocal(),
                                APPLY_AGS_SCATTER_SOURCES | APPLY_WGS_SCATTER_SOURCES |
                                  APPLY_AGS_FISSION_SOURCES | APPLY_WGS_FISSION_SOURCES);
  }

  return source_moments;
}

//###################################################################
/**Writes a given flux-moments vector to file.*/
void
lbs::LBSSolver::WriteFluxMoments(const std::string& file_base,
                                 const std::vector<double>& flux_moments)
{
  std::string file_name = file_base + std::to_string(Chi::mpi.location_id) + ".data";

  //============================================= Open file
  Chi::log.Log() << "Writing flux-moments to files with base-name " << file_base
                 << " and extension .data";
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
  std::string header_info = "Chi-Tech LinearBoltzmann: Flux moments file\n"
                            "Header size: 500 bytes\n"
                            "Structure(type-info):\n"
                            "uint64_t num_local_nodes\n"
                            "uint64_t num_moments\n"
                            "uint64_t num_groups\n"
                            "uint64_t num_records\n"
                            "uint64_t num_cells\n"
                            "Each cell:\n"
                            "  uint64_t cell_global_id\n"
                            "  uint64_t num_nodes\n"
                            "  Each node:\n"
                            "    double   x_position\n"
                            "    double   y_position\n"
                            "    double   z_position\n"
                            "Each record:\n"
                            "  uint64_t     cell_global_id\n"
                            "  unsigned int node_number\n"
                            "  unsigned int moment_num\n"
                            "  unsigned int group_num\n"
                            "  double       flux_moment_value\n";

  int header_size = (int)header_info.length();

  char header_bytes[500];
  memset(header_bytes, '-', 500);
  strncpy(header_bytes, header_info.c_str(), std::min(header_size, 499));
  header_bytes[499] = '\0';

  file << header_bytes;

  //============================================= Get relevant items
  auto NODES_ONLY = chi_math::UnknownManager::GetUnitaryUnknownManager();
  auto& sdm = discretization_;
  uint64_t num_local_nodes = discretization_->GetNumLocalDOFs(NODES_ONLY);
  uint64_t num_moments_t = static_cast<uint64_t>(num_moments_);
  uint64_t num_groups_t = static_cast<uint64_t>(num_groups_);
  uint64_t num_local_dofs = discretization_->GetNumLocalDOFs(flux_moments_uk_man_);
  uint64_t num_local_cells = grid_ptr_->local_cells.size();

  //============================================= Write num_ quantities
  file.write((char*)&num_local_nodes, sizeof(uint64_t));
  file.write((char*)&num_moments_t, sizeof(uint64_t));
  file.write((char*)&num_groups_t, sizeof(uint64_t));
  file.write((char*)&num_local_dofs, sizeof(uint64_t));
  file.write((char*)&num_local_cells, sizeof(uint64_t));

  //============================================= Write nodal positions for
  //                                              each cell
  for (const auto& cell : grid_ptr_->local_cells)
  {
    uint64_t cell_global_id = static_cast<uint64_t>(cell.global_id_);
    file.write((char*)&cell_global_id, sizeof(uint64_t));

    uint64_t num_nodes = discretization_->GetCellNumNodes(cell);
    file.write((char*)&num_nodes, sizeof(uint64_t));

    auto node_locations = discretization_->GetCellNodeLocations(cell);
    for (const auto& node : node_locations)
    {
      file.write((char*)&node.x, sizeof(double));
      file.write((char*)&node.y, sizeof(double));
      file.write((char*)&node.z, sizeof(double));
    } // for node
  }   // for cell

  //============================================= Write per dof data
  for (const auto& cell : grid_ptr_->local_cells)
    for (unsigned int i = 0; i < sdm->GetCellNumNodes(cell); ++i)
      for (unsigned int m = 0; m < num_moments_t; ++m)
        for (unsigned int g = 0; g < num_groups_; ++g)
        {
          uint64_t cell_global_id = cell.global_id_;
          uint64_t dof_map = sdm->MapDOFLocal(cell, i, flux_moments_uk_man_, m, g);

          assert(dof_map < flux_moments.size());
          double value = flux_moments[dof_map];

          file.write((char*)&cell_global_id, sizeof(uint64_t));
          file.write((char*)&i, sizeof(unsigned int));
          file.write((char*)&m, sizeof(unsigned int));
          file.write((char*)&g, sizeof(unsigned int));
          file.write((char*)&value, sizeof(double));
        }

  //============================================= Clean-up
  file.close();
}

//###################################################################
/**Reads a flux-moments vector from a file in the specified vector.*/
void
lbs::LBSSolver::ReadFluxMoments(const std::string& file_base,
                                std::vector<double>& flux_moments,
                                bool single_file /*=false*/)
{
  std::string file_name = file_base + std::to_string(Chi::mpi.location_id) + ".data";
  if (single_file) file_name = file_base + ".data";

  //============================================= Open file
  Chi::log.Log() << "Reading flux-moments file " << file_name;
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
  auto& sdm = discretization_;
  uint64_t num_local_nodes = discretization_->GetNumLocalDOFs(NODES_ONLY);
  uint64_t num_moments_t = static_cast<uint64_t>(num_moments_);
  uint64_t num_groups_t = static_cast<uint64_t>(num_groups_);
  uint64_t num_local_dofs = discretization_->GetNumLocalDOFs(flux_moments_uk_man_);
  uint64_t num_local_cells = grid_ptr_->local_cells.size();

  uint64_t file_num_local_nodes;
  uint64_t file_num_moments;
  uint64_t file_num_groups;
  uint64_t file_num_local_dofs;
  uint64_t file_num_local_cells;

  flux_moments.assign(num_local_dofs, 0.0);

  //============================================= Read header
  char header_bytes[500];
  header_bytes[499] = '\0';
  file.read(header_bytes, 499);

  file.read((char*)&file_num_local_nodes, sizeof(uint64_t));
  file.read((char*)&file_num_moments, sizeof(uint64_t));
  file.read((char*)&file_num_groups, sizeof(uint64_t));
  file.read((char*)&file_num_local_dofs, sizeof(uint64_t));
  file.read((char*)&file_num_local_cells, sizeof(uint64_t));

  //============================================= Check compatibility
  if (not single_file)
    if (file_num_local_nodes != num_local_nodes or file_num_moments != num_moments_t or
        file_num_groups != num_groups_t or file_num_local_dofs != num_local_dofs or
        file_num_local_cells != num_local_cells)
    {
      std::stringstream outstr;
      outstr << "num_local_nodes: " << file_num_local_nodes << " vs " << num_local_nodes << "\n";
      outstr << "num_moments_    : " << file_num_moments << " vs " << num_moments_t << "\n";
      outstr << "num_groups     : " << file_num_groups << " vs " << num_groups_ << "\n";
      outstr << "num_local_dofs : " << file_num_local_dofs << " vs " << num_local_dofs << "\n";
      outstr << "num_local_cells: " << file_num_local_cells << " vs " << num_local_cells << "\n";
      Chi::log.LogAll() << "Incompatible DOF data found in file " << file_name << "\n"
                        << "File data vs system:\n"
                        << outstr.str();
      file.close();
      return;
    }

  //============================================= Read cell nodal locations
  std::map<uint64_t, std::map<uint64_t, uint64_t>> file_cell_nodal_mapping;
  for (uint64_t c = 0; c < file_num_local_cells; ++c)
  {
    //============================ Read cell-id and num_nodes
    uint64_t cell_global_id;
    uint64_t num_nodes;

    file.read((char*)&cell_global_id, sizeof(uint64_t));
    file.read((char*)&num_nodes, sizeof(uint64_t));

    //============================ Read node locations
    std::vector<chi_mesh::Vector3> file_node_locations;
    file_node_locations.reserve(num_nodes);
    for (uint64_t n = 0; n < num_nodes; ++n)
    {
      double x, y, z;
      file.read((char*)&x, sizeof(double));
      file.read((char*)&y, sizeof(double));
      file.read((char*)&z, sizeof(double));

      file_node_locations.emplace_back(x, y, z);
    } // for file node n

    if (not grid_ptr_->IsCellLocal(cell_global_id)) continue;

    const auto& cell = grid_ptr_->cells[cell_global_id];

    //================ Now map file nodes to system nodes
    auto system_node_locations = discretization_->GetCellNodeLocations(cell);
    std::map<uint64_t, uint64_t> mapping;

    // Check num_nodes equal
    if (system_node_locations.size() != num_nodes)
      throw std::logic_error(std::string(__FUNCTION__) +
                             ": Incompatible number of nodes for a cell was encountered. Mapping "
                             "could not be performed.");

    bool mapping_successful = true; // Assume true, now try to disprove

    const auto& sys_nodes = system_node_locations;
    const auto& file_nodes = file_node_locations;
    size_t num_system_nodes = system_node_locations.size();

    for (uint64_t n = 0; n < num_nodes; ++n)
    {
      bool mapping_found = false;
      for (uint64_t m = 0; m < num_system_nodes; ++m)
        if ((sys_nodes[m] - file_nodes[n]).NormSquare() < 1.0e-12)
        {
          mapping[n] = m;
          mapping_found = true;
        }

      if (not mapping_found)
      {
        mapping_successful = false;
        break;
      }
    } // for n

    if (not mapping_successful)
      throw std::logic_error(std::string(__FUNCTION__) +
                             ": Incompatible node locations for a cell was encountered. Mapping "
                             "unsuccessful.");

    file_cell_nodal_mapping[cell_global_id] = std::move(mapping);
  } // for c (cell in file)

  //============================================= Commit to reading the file
  for (size_t dof = 0; dof < file_num_local_dofs; ++dof)
  {
    uint64_t cell_global_id;
    unsigned int node;
    unsigned int moment;
    unsigned int group;
    double flux_value;

    file.read((char*)&cell_global_id, sizeof(uint64_t));
    file.read((char*)&node, sizeof(unsigned int));
    file.read((char*)&moment, sizeof(unsigned int));
    file.read((char*)&group, sizeof(unsigned int));
    file.read((char*)&flux_value, sizeof(double));

    if (grid_ptr_->IsCellLocal(cell_global_id))
    {
      const auto& cell = grid_ptr_->cells[cell_global_id];
      const auto& node_mapping = file_cell_nodal_mapping.at(cell_global_id);

      size_t node_mapped = node_mapping.at(node);

      size_t dof_map = sdm->MapDOFLocal(cell, node_mapped, flux_moments_uk_man_, moment, group);

      assert(dof_map < flux_moments.size());
      flux_moments[dof_map] = flux_value;
    } // if cell is local
  }   // for dof

  //============================================= Clean-up
  file.close();
}