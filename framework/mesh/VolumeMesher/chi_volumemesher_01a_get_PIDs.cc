#include "mesh/VolumeMesher/chi_volumemesher.h"

#include "mesh/MeshHandler/chi_meshhandler.h"
#include "mesh/VolumeMesher/Extruder/volmesher_extruder.h"
#include "mesh/VolumeMesher/PredefinedUnpartitioned/volmesher_predefunpart.h"

#include "chi_runtime.h"
#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/**Obtains the xy partition IDs of a cell.
 * Cell xy_partition ids are obtained from
 * the surface mesher.*/
std::pair<int, int>
chi_mesh::VolumeMesher::GetCellXYPartitionID(chi_mesh::Cell* cell)
{
  std::pair<int, int> ij_id(0, 0);

  if (Chi::mpi.process_count == 1) { return ij_id; }

  //================================================== Get the current handler
  auto& mesh_handler = chi_mesh::GetCurrentHandler();
  auto& vol_mesher = mesh_handler.GetVolumeMesher();

  //====================================== Sanity check on partitioning
  size_t num_x_subsets = vol_mesher.options.xcuts.size() + 1;
  size_t num_y_subsets = vol_mesher.options.ycuts.size() + 1;

  size_t x_remainder = num_x_subsets % vol_mesher.options.partition_x;
  size_t y_remainder = num_y_subsets % vol_mesher.options.partition_y;

  if (x_remainder != 0)
  {
    Chi::log.LogAllError() << "When specifying x-partitioning, the number of grp_subsets in x "
                              "needs to be divisible by the number of partitions in x.";
    Chi::Exit(EXIT_FAILURE);
  }

  if (y_remainder != 0)
  {
    Chi::log.LogAllError() << "When specifying y-partitioning, the number of grp_subsets in y "
                              "needs to be divisible by the number of partitions in y.";
    Chi::Exit(EXIT_FAILURE);
  }

  size_t subsets_per_partitionx = num_x_subsets / vol_mesher.options.partition_x;
  size_t subsets_per_partitiony = num_y_subsets / vol_mesher.options.partition_y;

  //====================================== Determine x-partition
  int x = -1;
  int xcount = -1;
  for (size_t i = subsets_per_partitionx - 1; i < vol_mesher.options.xcuts.size();
       i += subsets_per_partitionx)
  {
    xcount++;
    if (cell->centroid_.x <= vol_mesher.options.xcuts[i])
    {
      x = xcount;
      break;
    }
  }
  if (x < 0) { x = vol_mesher.options.partition_x - 1; }

  //====================================== Determine y-partition
  int y = -1;
  int ycount = -1;
  for (size_t i = subsets_per_partitiony - 1; i < vol_mesher.options.ycuts.size();
       i += subsets_per_partitiony)
  {
    ycount++;
    if (cell->centroid_.y <= vol_mesher.options.ycuts[i])
    {
      y = ycount;
      break;
    }
  }
  if (y < 0) { y = vol_mesher.options.partition_y - 1; }

  //====================================== Set partitioning
  ij_id.first = x;
  ij_id.second = y;

  return ij_id;
}

//###################################################################
/**Obtains the xyz partition IDs of a cell.
 * Cell xy_partition ids are obtained from
 * the surface mesher. z id is obtained from the volume mesher.*/
std::tuple<int, int, int>
chi_mesh::VolumeMesher::GetCellXYZPartitionID(chi_mesh::Cell* cell)
{
  std::tuple<int, int, int> ijk_id(0, 0, 0);
  bool found_partition = false;

  if (Chi::mpi.process_count == 1) { return ijk_id; }

  //================================================== Get ij indices
  std::pair<int, int> ij_id = GetCellXYPartitionID(cell);

  //================================================== Get the current handler
  auto& mesh_handler = chi_mesh::GetCurrentHandler();
  auto& vol_mesher = mesh_handler.GetVolumeMesher();

  if (vol_mesher.options.partition_z == 1)
  {
    found_partition = true;
    std::get<0>(ijk_id) = ij_id.first;
    std::get<1>(ijk_id) = ij_id.second;
    std::get<2>(ijk_id) = 0;
  }
  else if (vol_mesher.Type() == VolumeMesherType::EXTRUDER)
  {
    auto extruder = dynamic_cast<chi_mesh::VolumeMesherExtruder&>(vol_mesher);
    const auto& vertex_layers = extruder.GetVertexLayers();
    //====================================== Create virtual cuts
    if (vol_mesher.options.zcuts.empty())
    {
      size_t num_sub_layers = vertex_layers.size() - 1;

      if ((num_sub_layers % vol_mesher.options.partition_z) != 0)
      {
        Chi::log.LogAllError() << "Number of sub-layers in extruded mesh is not divisible "
                               << "by the requested number of z-partitions.";
        Chi::Exit(EXIT_FAILURE);
      }

      int delta_zk = num_sub_layers / vol_mesher.options.partition_z;
      for (int k = 0; k < (vol_mesher.options.partition_z); k++)
      {
        int layer_index = k * delta_zk + delta_zk;
        if (layer_index > (vertex_layers.size() - 1))
        {
          layer_index = (int)vertex_layers.size() - 1;
          vol_mesher.options.zcuts.push_back(vertex_layers[layer_index]);
        }
        else
        {
          vol_mesher.options.zcuts.push_back(vertex_layers[layer_index]);

          if (Chi::log.GetVerbosity() == chi::ChiLog::LOG_LVL::LOG_0VERBOSE_2)
          {
            printf("Z-Cut %lu, %g\n", vol_mesher.options.zcuts.size(), vertex_layers[layer_index]);
          }
        }
      }
    }

    //====================================== Scan cuts for location
    double zmin = -1.0e-16;
    for (int k = 0; k < (vol_mesher.options.zcuts.size()); k++)
    {
      double zmax = vol_mesher.options.zcuts[k];

      double z = cell->centroid_.z;

      if (Chi::log.GetVerbosity() == chi::ChiLog::LOG_0VERBOSE_2)
      {
        printf("zmax = %g, zmin = %g, cell_z = %g\n", zmax, zmin, z);
      }

      if ((z > zmin) && (z < zmax))
      {
        std::get<0>(ijk_id) = ij_id.first;
        std::get<1>(ijk_id) = ij_id.second;
        std::get<2>(ijk_id) = k;

        found_partition = true;
        break;
      }
      zmin = zmax;
    }
  } // if typeid
  else if (vol_mesher.Type() == VolumeMesherType::UNPARTITIONED)
  {
    if (vol_mesher.options.zcuts.empty())
    {
      throw std::invalid_argument("Cell z-partitioning cannot be determined "
                                  "because no z-cuts are supplied to volume "
                                  "mesher.");
    }

    //====================================== Scan cuts for location
    std::vector<double> temp_zcuts = vol_mesher.options.zcuts;
    double zmin = -1.0e16;
    double zmax = 1.0e16;
    temp_zcuts.push_back(zmax);
    for (int k = 0; k < (temp_zcuts.size()); k++)
    {
      zmax = temp_zcuts[k];

      double z = cell->centroid_.z;

      if (Chi::log.GetVerbosity() == chi::ChiLog::LOG_0VERBOSE_2)
      {
        printf("zmax = %g, zmin = %g, cell_z = %g\n", zmax, zmin, z);
      }

      if ((z > zmin) && (z < zmax))
      {
        std::get<0>(ijk_id) = ij_id.first;
        std::get<1>(ijk_id) = ij_id.second;
        std::get<2>(ijk_id) = k;

        found_partition = true;
        break;
      }
      zmin = zmax;
    } // for k
  }   // if typeid

  //================================================== Report unallocated
  // item_id
  if (!found_partition)
  {
    Chi::log.LogAllError() << "A cell was encountered for which "
                              "no zpartition id was found";
    Chi::Exit(EXIT_FAILURE);
  }

  return ijk_id;
}