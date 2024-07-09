// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/materials/isotropic_multigroup_source.h"
#include "framework/math/math.h"
#include <functional>
#include <map>

namespace opensn
{

using DirIDs = std::vector<size_t>; ///< Direction-IDs
using UniqueSOGroupings = std::vector<DirIDs>;
using DirIDToSOMap = std::map<size_t, size_t>;
using MatVec3 = std::vector<std::vector<Vector3>>;

enum class SolverType
{
  DISCRETE_ORDINATES = 1,
  DIFFUSION_DFEM = 2,
  DIFFUSION_CFEM = 3,
};

enum class GeometryType
{
  NO_GEOMETRY_SET = 0,
  ONED_SLAB = 1,
  ONED_CYLINDRICAL = 2,
  ONED_SPHERICAL = 3,
  TWOD_CARTESIAN = 4,
  TWOD_CYLINDRICAL = 5,
  THREED_CARTESIAN = 6
};

inline CoordinateSystemType
MapGeometryTypeToCoordSys(const GeometryType gtype)
{
  switch (gtype)
  {
    case GeometryType::ONED_SLAB:
    case GeometryType::TWOD_CARTESIAN:
    case GeometryType::THREED_CARTESIAN:
      return CoordinateSystemType::CARTESIAN;
    case GeometryType::ONED_SPHERICAL:
      return CoordinateSystemType::SPHERICAL;
    case GeometryType::ONED_CYLINDRICAL:
    case GeometryType::TWOD_CYLINDRICAL:
      return CoordinateSystemType::CYLINDRICAL;
    default:
      return CoordinateSystemType::CARTESIAN;
  }
}

enum class AngleAggregationType
{
  UNDEFINED = 0,
  SINGLE = 1,
  POLAR = 2,
  AZIMUTHAL = 3,
};

enum class BoundaryType
{
  VACUUM = 1,     ///< Zero for all angles, space
  ISOTROPIC = 2,  ///< One value for all angles, homogenous in space
  REFLECTING = 3, ///< Reflecting boundary condition about a normal
  ARBITRARY = 4   ///< Complex different for each angle and face node
};

struct BoundaryPreference
{
  BoundaryType type;
  std::vector<double> isotropic_mg_source;
  std::string source_function;
};

enum SourceType
{
  APPLY_FIXED_SOURCES = (1 << 0),
  APPLY_WGS_SCATTER_SOURCES = (1 << 1),
  APPLY_AGS_SCATTER_SOURCES = (1 << 2),
  APPLY_WGS_FISSION_SOURCES = (1 << 3),
  APPLY_AGS_FISSION_SOURCES = (1 << 4),
  SUPPRESS_WG_SCATTER = (1 << 5),
  ZERO_INCOMING_DELAYED_PSI = (1 << 6)
};

/**
 * SourceFlags is a combination of `SourceType`s
 */
struct SourceFlags
{
  /**
   * Create an empty `SourceFlags`
   */
  SourceFlags() : flags_(0) {}

  /**
   * Create `SourceFlags` from a `SourceType`
   *
   * \param type Type to set in the `SourceFlags`
   */
  SourceFlags(SourceType type) : flags_(type) {}

  /**
   * Add flags using the `|=` operator
   *
   * \param type Type to set
   * \return Combination of our types and new `type`
   */
  SourceFlags& operator|=(SourceType type)
  {
    flags_ |= type;
    return *this;
  }

  /**
   * Add source using the `|=` operator
   *
   * \param src Source to add
   * \return Combination of our types and new types from `src`
   */
  SourceFlags& operator|=(const SourceFlags& src)
  {
    flags_ |= src.flags_;
    return *this;
  }

  /**
   * Test is there are any types set
   * \return `true` if there are no types set, `false` otherwise
   */
  bool Empty() const { return flags_ == 0; }

  /**
   * Unset a type
   *
   * \param type Type to unset
   */
  void Unset(SourceType type) { flags_ &= ~type; }

  /**
   * Test if type is set using the `&` operator
   *
   * \param type Type to test
   * \return `true` if `type` is set, `false` otherwise
   */
  bool operator&(const SourceType& type) const { return flags_ & type; }

private:
  /// Combination of `SourceFlags`
  int flags_;
};

/**
 * Operator to join individual `SourceFlags`s using the `|` operator
 *
 * \param s1 First source
 * \param s2 Second source
 * \return Source with types from both sources
 */
inline SourceFlags
operator|(SourceFlags s1, SourceFlags s2)
{
  SourceFlags src = s1;
  src |= s2;
  return src;
}

/**
 * Two `SourceType`s get combined into `SourceFlags`
 *
 * \param f1 First source type
 * \param f2 Second source type
 * \return Source with both types set
 */
inline SourceFlags
operator|(SourceType f1, SourceType f2)
{
  SourceFlags src;
  src |= f1;
  src |= f2;
  return src;
}

enum class PhiSTLOption
{
  PHI_OLD = 1,
  PHI_NEW = 2
};

class LBSGroupset;

typedef std::function<void(const LBSGroupset& groupset,
                           std::vector<double>& q,
                           const std::vector<double>& phi,
                           const std::vector<double>& densities,
                           const SourceFlags source_flags)>
  SetSourceFunction;

/**Struct for storing LBS options.*/
struct LBSOptions
{
  GeometryType geometry_type = GeometryType::NO_GEOMETRY_SET;
  SpatialDiscretizationType sd_type = SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS;
  unsigned int scattering_order = 1;
  int max_mpi_message_size = 32768;

  std::filesystem::path read_restart_path;
  std::filesystem::path write_restart_path =
    opensn::input_path.replace_extension("restart").string() + "/" +
    opensn::input_path.stem().string();
  size_t write_restart_time_interval = 0;

  bool enable_ags_restart_write = true;

  bool use_precursors = false;
  bool use_src_moments = false;

  bool save_angular_flux = false;

  bool adjoint = false;

  bool verbose_inner_iterations = true;
  bool verbose_ags_iterations = false;
  bool verbose_outer_iterations = true;
  int max_ags_iterations = 1;
  double ags_tolerance = 1.0e-6;

  bool power_field_function_on = false;
  double power_default_kappa = 3.20435e-11; // 200MeV to Joule
  double power_normalization = -1.0;

  std::string field_function_prefix_option = "prefix";
  std::string field_function_prefix; // Default is empty

  LBSOptions() = default;
};

/**Transport view of a cell*/
class CellLBSView
{
private:
  size_t phi_address_;
  int num_nodes_;
  int num_groups_;
  int num_grps_moms_;
  const MultiGroupXS* xs_;
  double volume_;
  const std::vector<bool> face_local_flags_;
  const std::vector<int> face_locality_;
  const std::vector<const Cell*> neighbor_cell_ptrs_;
  std::vector<std::vector<double>> outflow_;

public:
  CellLBSView(size_t phi_address,
              int num_nodes,
              int num_groups,
              int num_moments,
              int num_faces,
              const MultiGroupXS& xs_mapping,
              double volume,
              const std::vector<bool>& face_local_flags,
              const std::vector<int>& face_locality,
              const std::vector<const Cell*>& neighbor_cell_ptrs,
              bool cell_on_boundary)
    : phi_address_(phi_address),
      num_nodes_(num_nodes),
      num_groups_(num_groups),
      num_grps_moms_(num_groups * num_moments),
      xs_(&xs_mapping),
      volume_(volume),
      face_local_flags_(face_local_flags),
      face_locality_(face_locality),
      neighbor_cell_ptrs_(neighbor_cell_ptrs)
  {
    if (cell_on_boundary)
      outflow_.resize(num_faces);
  }

  size_t MapDOF(int node, int moment, int grp) const
  {
    return phi_address_ + node * num_grps_moms_ + num_groups_ * moment + grp;
  }

  const MultiGroupXS& XS() const { return *xs_; }

  bool IsFaceLocal(int f) const { return face_local_flags_[f]; }

  int FaceLocality(int f) const { return face_locality_[f]; }

  const Cell* FaceNeighbor(int f) const { return neighbor_cell_ptrs_[f]; }

  int NumNodes() const { return num_nodes_; }

  double Volume() const { return volume_; }

  void ZeroOutflow(int f, int g)
  {
    if (f < outflow_.size() and g < outflow_[f].size())
      outflow_[f][g] = 0.0;
  }

  void AddOutflow(int f, int g, double intS_mu_psi)
  {
    if (f < outflow_.size())
    {
      if (outflow_[f].empty())
        outflow_[f].resize(num_groups_, 0.0);
      if (g < outflow_[f].size())
        outflow_[f][g] += intS_mu_psi;
    }
  }

  double GetOutflow(int f, int g) const
  {
    if (f < outflow_.size() and g < outflow_[f].size())
      return outflow_[f][g];
    return 0.0;
  }

  void ReassignXS(const MultiGroupXS& xs) { xs_ = &xs; }
};

struct UnitCellMatrices
{
  std::vector<std::vector<double>> intV_gradshapeI_gradshapeJ;
  std::vector<std::vector<Vector3>> intV_shapeI_gradshapeJ;
  std::vector<std::vector<double>> intV_shapeI_shapeJ;
  std::vector<double> intV_shapeI;

  std::vector<std::vector<std::vector<double>>> intS_shapeI_shapeJ;
  std::vector<std::vector<std::vector<Vector3>>> intS_shapeI_gradshapeJ;
  std::vector<std::vector<double>> intS_shapeI;
};

} // namespace opensn
