#pragma once

#include "B_DiscreteOrdinatesSolver/lbs_discrete_ordinates_solver.h"
#include "LinearBoltzmannSolvers/A_LBSSolver/Groupset/lbs_groupset.h"

namespace lbs
{

/** A neutral particle transport solver in point-symmetric and axial-symmetric
 *  curvilinear coordinates. */
class DiscreteOrdinatesCurvilinearSolver : public lbs::DiscreteOrdinatesSolver
{
  //  Attributes
private:
  /** Coordinate system type. */
  chi_math::CoordinateSystemType coord_system_type_;
  /** Discretisation pointer to matrices of the secondary cell view
   *  (matrices of the primary cell view forwarded to the base class). */
  std::shared_ptr<chi_math::SpatialDiscretization> discretization_secondary_;
  std::vector<lbs::UnitCellMatrices> secondary_unit_cell_matrices_;

  //  Methods
public:
  static chi::InputParameters GetInputParameters();
  explicit DiscreteOrdinatesCurvilinearSolver(const chi::InputParameters& params);

  // Deleted copy constructors
  DiscreteOrdinatesCurvilinearSolver(const DiscreteOrdinatesCurvilinearSolver&) = delete;
  DiscreteOrdinatesCurvilinearSolver& operator=(const DiscreteOrdinatesCurvilinearSolver&) = delete;

protected:
  void PerformInputChecks() override;
  void InitializeSpatialDiscretization() override;
  void ComputeSecondaryUnitIntegrals();

private:
  std::shared_ptr<SweepChunk> SetSweepChunk(lbs::LBSGroupset& groupset) override;
};

} // namespace lbs
