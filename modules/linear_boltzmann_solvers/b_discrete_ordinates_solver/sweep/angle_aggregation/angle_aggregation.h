#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/sweep_namespace.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/spds/spds.h"
#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/angle_set/angle_set_group.h"
#include "framework/math/quadratures/angular_quadrature_base.h"
#include <memory>

namespace opensn
{
namespace lbs
{

/**Angle aggregation has to cater for running the 8 corners of a 3D
 * partitioning, the 4 corners of a 2D partitioning (the latter 2 both being
 * polar angle aggregation) as well as single angle aggregation.
 *
 * At the most fundamental level this manifests as a number of angle indices
 * that share a SPDS, however SPDS do not have to be unique which allows
 * the notion of polar angle sets.
 * For single angle aggregation each SPDS is associated
 * with a single angle index. For polar angle aggregation a single SPDS can
 * have multiple indices associated with the same azimuthal angle but
 * different polar angles. We call this manifestation a "AngleSet".
 *
 * The octant based separation is achieved via the notion of "AngleSetGroup"
 * which will group angle sets for each quadrant or octant
 * (depending on 2D or 3D).*/
class AngleAggregation
{
public:
  std::vector<AngleSetGroup> angle_set_groups;
  std::map<uint64_t, std::shared_ptr<SweepBndry>> sim_boundaries;
  size_t number_of_groups = 0;
  size_t number_of_group_subsets = 0;
  std::shared_ptr<AngularQuadrature> quadrature = nullptr;

private:
  bool is_setup = false;
  std::pair<size_t, size_t> number_angular_unknowns;
  bool num_ang_unknowns_avail = false;

public:
  std::shared_ptr<MeshContinuum> grid = nullptr;

  /** Sets up the angle-aggregation object. */
  AngleAggregation(const std::map<uint64_t, std::shared_ptr<SweepBndry>>& in_sim_boundaries,
                   size_t in_number_of_groups,
                   size_t in_number_of_group_subsets,
                   std::shared_ptr<AngularQuadrature>& in_quadrature,
                   std::shared_ptr<MeshContinuum>& in_grid);

  bool IsSetup() const { return is_setup; }

public:
  /** Resets all the outgoing intra-location and inter-location
   * cyclic interfaces.*/
  void ZeroOutgoingDelayedPsi();
  /** Resets all the incoming intra-location and inter-location
   * cyclic interfaces.*/
  void ZeroIncomingDelayedPsi();

  /** Initializes reflecting boundary conditions. */
  void InitializeReflectingBCs();

  /** Returns a pair of numbers containing the number of
   * delayed angular unknowns both locally and globally, respectively. */
  std::pair<size_t, size_t> GetNumDelayedAngularDOFs();

  /** Assembles angular unknowns into the reference vector. */
  void AppendNewDelayedAngularDOFsToArray(int64_t& index, double* x_ref);
  /** Assembles angular unknowns into the reference vector. */
  void AppendOldDelayedAngularDOFsToArray(int64_t& index, double* x_ref);

  /** Assembles angular unknowns into the reference vector. */
  void SetOldDelayedAngularDOFsFromArray(int64_t& index, const double* x_ref);
  /** Assembles angular unknowns into the reference vector. */
  void SetNewDelayedAngularDOFsFromArray(int64_t& index, const double* x_ref);

  /**Gets the current values of the angular unknowns as an STL vector.*/
  std::vector<double> GetNewDelayedAngularDOFsAsSTLVector();
  /**Gets the current values of the angular unknowns as an STL vector.*/
  void SetNewDelayedAngularDOFsFromSTLVector(const std::vector<double>& stl_vector);

  /**Gets the current values of the angular unknowns as an STL vector.*/
  std::vector<double> GetOldDelayedAngularDOFsAsSTLVector();
  /**Gets the current values of the angular unknowns as an STL vector.*/
  void SetOldDelayedAngularDOFsFromSTLVector(const std::vector<double>& stl_vector);

  /**Copies the old delayed angular fluxes to the new.*/
  void SetDelayedPsiOld2New();
  /**Copies the new delayed angular fluxes to the old.*/
  void SetDelayedPsiNew2Old();
};

} // namespace lbs
} // namespace opensn
