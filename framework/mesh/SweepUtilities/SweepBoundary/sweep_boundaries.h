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

//###################################################################
/**Base class for sweep related boundaries.*/
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

  virtual double* HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                           unsigned int face_num,
                                           unsigned int fi,
                                           unsigned int angle_num,
                                           int group_num,
                                           size_t gs_ss_begin);

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

//###################################################################
/** Zero fluxes homogenous on a boundary and in angle.*/
class BoundaryVaccuum : public SweepBoundary
{
private:
  std::vector<double> boundary_flux_;

public:
  explicit BoundaryVaccuum(
    size_t in_num_groups,
    chi_math::CoordinateSystemType coord_type = chi_math::CoordinateSystemType::CARTESIAN)
    : SweepBoundary(BoundaryType::INCIDENT_VACCUUM, in_num_groups, coord_type),
      boundary_flux_(in_num_groups, 0.0)
  {
  }

  double* HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                   unsigned int face_num,
                                   unsigned int fi,
                                   unsigned int angle_num,
                                   int group_num,
                                   size_t gs_ss_begin) override;
};

//###################################################################
/** Specified incident fluxes homogenous on a boundary.*/
class BoundaryIsotropicHomogenous : public SweepBoundary
{
private:
  std::vector<double> boundary_flux;

public:
  explicit BoundaryIsotropicHomogenous(
    size_t in_num_groups,
    std::vector<double> ref_boundary_flux,
    chi_math::CoordinateSystemType coord_type = chi_math::CoordinateSystemType::CARTESIAN)
    : SweepBoundary(BoundaryType::INCIDENT_ISOTROPIC_HOMOGENOUS, in_num_groups, coord_type),
      boundary_flux(std::move(ref_boundary_flux))
  {
  }

  double* HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                   unsigned int face_num,
                                   unsigned int fi,
                                   unsigned int angle_num,
                                   int group_num,
                                   size_t gs_ss_begin) override;
};

//###################################################################
/** Reflective boundary condition.*/
class BoundaryReflecting : public SweepBoundary
{
protected:
  const chi_mesh::Normal normal_;
  bool opposing_reflected_ = false;

  typedef std::vector<double> DOFVec;   // Groups per DOF
  typedef std::vector<DOFVec> FaceVec;  // DOFs per face
  typedef std::vector<FaceVec> CellVec; // Faces per cell
  typedef std::vector<CellVec> AngVec;  // Cell per angle

  // angle,cell,face,dof,group
  // Populated by angle aggregation
  std::vector<AngVec> hetero_boundary_flux_;
  std::vector<AngVec> hetero_boundary_flux_old_;

  std::vector<int> reflected_anglenum_;
  std::vector<std::vector<bool>> angle_readyflags_;

public:
  BoundaryReflecting(
    size_t in_num_groups,
    const chi_mesh::Normal& in_normal,
    chi_math::CoordinateSystemType coord_type = chi_math::CoordinateSystemType::CARTESIAN)
    : SweepBoundary(BoundaryType::REFLECTING, in_num_groups, coord_type), normal_(in_normal)
  {
  }

  const chi_mesh::Vector3& Normal() const { return normal_; }
  bool IsOpposingReflected() const { return opposing_reflected_; }
  void SetOpposingReflected(bool value) { opposing_reflected_ = value; }

  std::vector<AngVec>& GetHeteroBoundaryFluxNew() { return hetero_boundary_flux_; }
  std::vector<AngVec>& GetHeteroBoundaryFluxOld() { return hetero_boundary_flux_old_; }

  std::vector<int>& GetReflectedAngleIndexMap() { return reflected_anglenum_; }
  std::vector<std::vector<bool>>& GetAngleReadyFlags() { return angle_readyflags_; }

  double* HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                   unsigned int face_num,
                                   unsigned int fi,
                                   unsigned int angle_num,
                                   int group_num,
                                   size_t gs_ss_begin) override;
  double* HeterogeneousPsiOutgoing(uint64_t cell_local_id,
                                   unsigned int face_num,
                                   unsigned int fi,
                                   unsigned int angle_num,
                                   size_t gs_ss_begin) override;

  void UpdateAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss) override;
  bool CheckAnglesReadyStatus(const std::vector<size_t>& angles, size_t gs_ss) override;
  void ResetAnglesReadyStatus();
};

/**This boundary function class can be derived from to
 * provide a much more custom experience. This function
 * is called during Setup. */
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

//###################################################################
/** Specified incident fluxes homogenous on a boundary.*/
class BoundaryIncidentHeterogeneous : public SweepBoundary
{
private:
  std::unique_ptr<BoundaryFunction> boundary_function_;
  const uint64_t ref_boundary_id_;

  typedef std::vector<double> FaceNodeData;
  typedef std::vector<FaceNodeData> FaceData;
  typedef std::vector<FaceData> CellData;

  std::vector<CellData> local_cell_data_;

public:
  explicit BoundaryIncidentHeterogeneous(
    size_t in_num_groups,
    std::unique_ptr<BoundaryFunction> in_bndry_function,
    uint64_t in_ref_boundary_id,
    chi_math::CoordinateSystemType coord_type = chi_math::CoordinateSystemType::CARTESIAN)
    : SweepBoundary(BoundaryType::INCIDENT_ANISOTROPIC_HETEROGENEOUS, in_num_groups, coord_type),
      boundary_function_(std::move(in_bndry_function)),
      ref_boundary_id_(in_ref_boundary_id)
  {
  }

  double* HeterogeneousPsiIncoming(uint64_t cell_local_id,
                                   unsigned int face_num,
                                   unsigned int fi,
                                   unsigned int angle_num,
                                   int group_num,
                                   size_t gs_ss_begin) override;

  void Setup(const chi_mesh::MeshContinuum& grid,
             const chi_math::AngularQuadrature& quadrature) override;
};

} // namespace chi_mesh::sweep_management
