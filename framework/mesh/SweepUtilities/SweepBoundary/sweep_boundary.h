#pragma once

#include "mesh/chi_mesh.h"
#include "math/chi_math.h"

#include <vector>
#include <limits>

namespace chi_mesh::sweep_management
{

enum class BoundaryType
{
  INCIDENT_VACCUUM = 0,                  ///< Zero for all angles, space
  INCIDENT_ISOTROPIC_HOMOGENOUS = 1,     ///< One value for all angles, homogenous in space
  REFLECTING = 2,                        ///< Reflecting boundary condition about a normal
  INCIDENT_ANISOTROPIC_HETEROGENEOUS = 3 ///< Complex different for each angle and face node
};

/**
 * Base class for sweep related boundaries.
 */
class SweepBoundary
{
private:
  const chi_mesh::sweep_management::BoundaryType type_;
  const chi_math::CoordinateSystemType coord_type_;
  double evaluation_time_ = 0.0; ///< Time value passed to boundary functions
protected:
  std::vector<double> zero_boundary_flux_;
  size_t num_groups_;

public:
  explicit SweepBoundary(BoundaryType bndry_type,
                         size_t in_num_groups,
                         chi_math::CoordinateSystemType coord_type)
    : type_(bndry_type), coord_type_(coord_type), num_groups_(in_num_groups)
  {
    zero_boundary_flux_.resize(num_groups_, 0.0);
  }

  virtual ~SweepBoundary() = default;
  BoundaryType Type() const { return type_; }
  chi_math::CoordinateSystemType CoordType() const { return coord_type_; }
  bool IsReflecting() const { return type_ == BoundaryType::REFLECTING; }

  double GetEvaluationTime() const { return evaluation_time_; }
  void SetEvaluationTime(double time) { evaluation_time_ = time; }

  /**
   * Returns a pointer to a heterogeneous flux storage location.
   */
  virtual double* HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                           unsigned int face_num,
                                           unsigned int fi,
                                           unsigned int angle_num,
                                           int group_num,
                                           size_t gs_ss_begin);

  /**
   * Returns a pointer to a heterogeneous flux storage location.
   */
  virtual double* HeterogeneousPsiOutgoing(uint64_t cell_local_id,
                                           unsigned int face_num,
                                           unsigned int fi,
                                           unsigned int angle_num,
                                           size_t gs_ss_begin);

  virtual void UpdateAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss) {}

  virtual bool CheckAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss)
  {
    return true;
  }

  virtual void Setup(const chi_mesh::MeshContinuum& grid,
                     const chi_math::AngularQuadrature& quadrature)
  {
  }

  double* ZeroFlux(int group_num) { return &zero_boundary_flux_[group_num]; }
};

/**
 * This boundary function class can be derived from to
 * provide a much more custom experience. This function
 * is called during Setup.
 */
class BoundaryFunction
{
public:
  virtual std::vector<double>
  Evaluate(size_t cell_global_id,
           int cell_material_id,
           unsigned int face_index,
           unsigned int face_node_index,
           const chi_mesh::Vector3& face_node_location,
           const chi_mesh::Vector3& face_node_normal,
           const std::vector<int>& quadrature_angle_indices,
           const std::vector<chi_mesh::Vector3>& quadrature_angle_vectors,
           const std::vector<std::pair<double, double>>& quadrature_phi_theta_angles,
           const std::vector<int>& group_indices,
           double time) = 0;

  virtual ~BoundaryFunction() = default;
};

} // namespace chi_mesh::sweep_management
