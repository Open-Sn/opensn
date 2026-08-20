// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "hdf5.h"
#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace opensn
{

class DiscreteOrdinatesProblem;

struct UncollidedFluxData
{
  unsigned int num_groups = 0;
  unsigned int max_moment_order = 0;
  std::uint64_t global_cell_count = 0;
  double source_rate = 0.0;
  double outflow_rate = 0.0;
  std::vector<double> local_flux_moments;
};

class DiscreteOrdinatesProblemIO
{
public:
  static UncollidedFluxData ReadUncollidedFlux(const DiscreteOrdinatesProblem& do_problem,
                                               const std::string& file_name);

  static bool ReadRestartData(DiscreteOrdinatesProblem& do_problem,
                              hid_t file_id,
                              bool allow_transient_initialization_from_steady);

  static bool WriteRestartData(const DiscreteOrdinatesProblem& do_problem, hid_t file_id);

  /**
   * Write an angular flux vector to a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param opt_src Optional angular flux source vector
   */
  static void WriteAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    std::optional<const std::reference_wrapper<std::vector<std::vector<double>>>> opt_src =
      std::nullopt);

  /**
   * Read an angular flux vector from a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param opt_dest Optional angular flux destination vector
   */
  static void ReadAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    std::optional<std::reference_wrapper<std::vector<std::vector<double>>>> opt_dest =
      std::nullopt);

  /**
   * Surface Angular flux
   */
  using QuantizedCoordinate = std::tuple<int64_t, int64_t, int64_t>;
  using FaceMap = std::map<QuantizedCoordinate, uint64_t>;

  struct SurfaceMap
  {
    std::vector<uint64_t> cell_ids;
    std::vector<uint64_t> num_face_nodes;
    std::map<QuantizedCoordinate, uint64_t> cell_map;
    std::vector<uint64_t> cell_stride;
    std::vector<FaceMap> faces;
    std::vector<double> nodes_x;
    std::vector<double> nodes_y;
    std::vector<double> nodes_z;
  };

  struct SurfaceData
  {
    std::vector<double> omega;
    std::vector<double> mu;
    std::vector<double> wt_d;
    std::vector<double> mass_matrix;
    std::vector<double> fe_shape;
    std::vector<double> psi;
    std::vector<uint64_t> node_index;
    std::vector<uint64_t> dir_index;
  };

  struct SurfaceAngularFlux
  {
    int groupset_id = 0;
    std::string surface_name;
    SurfaceMap mapping;
    SurfaceData data;
  };

  /**
   * Write surface angular flux vector(s) to a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param bndry_surfs Boundary surface names
   * \param int_surfs Internal surface definitions
   */
  static void WriteSurfaceAngularFluxes(
    DiscreteOrdinatesProblem& do_problem,
    const std::string& file_base,
    const std::vector<std::string>& bndry_surfs,
    const std::vector<std::pair<std::string, std::pair<std::string, double>>>& int_surfs);

  /**
   * Read a surface angular flux vector from a file.
   *
   * \param do_problem Discrete ordinates problem
   * \param file_base File name base
   * \param surfaces Surface names to read
   */
  static std::vector<SurfaceAngularFlux>
  ReadSurfaceAngularFluxes(DiscreteOrdinatesProblem& do_problem,
                           const std::string& file_base,
                           const std::vector<std::string>& surfaces);
};

} // namespace opensn
