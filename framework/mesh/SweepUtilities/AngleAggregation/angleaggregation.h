#ifndef CHI_ANGLEAGGREGATION_H
#define CHI_ANGLEAGGREGATION_H

#include "mesh/SweepUtilities/sweep_namespace.h"
#include "mesh/SweepUtilities/SPDS/SPDS.h"
#include "mesh/SweepUtilities/AngleSetGroup/anglesetgroup.h"

#include "math/Quadratures/angular_quadrature_base.h"

#include <memory>

//###################################################################
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
class chi_mesh::sweep_management::AngleAggregation
{
public:
  typedef std::shared_ptr<SweepBndry> SweepBndryPtr;

public:
  std::vector<AngleSetGroup> angle_set_groups;
  std::map<uint64_t, SweepBndryPtr> sim_boundaries;
  size_t number_of_groups = 0;
  size_t number_of_group_subsets = 0;
  std::shared_ptr<chi_math::AngularQuadrature> quadrature = nullptr;

private:
  bool is_setup = false;
  std::pair<size_t, size_t> number_angular_unknowns;
  bool num_ang_unknowns_avail = false;

public:
  chi_mesh::MeshContinuumPtr grid = nullptr;

  AngleAggregation(const std::map<uint64_t, SweepBndryPtr>& in_sim_boundaries,
                   size_t in_number_of_groups,
                   size_t in_number_of_group_subsets,
                   std::shared_ptr<chi_math::AngularQuadrature>& in_quadrature,
                   chi_mesh::MeshContinuumPtr& in_grid);

  bool IsSetup() const { return is_setup; }

public:
  void ZeroOutgoingDelayedPsi();
  void ZeroIncomingDelayedPsi();

  void InitializeReflectingBCs();

  std::pair<size_t, size_t> GetNumDelayedAngularDOFs();

  void AppendNewDelayedAngularDOFsToArray(int64_t& index, double* x_ref);
  void AppendOldDelayedAngularDOFsToArray(int64_t& index, double* x_ref);

  void SetOldDelayedAngularDOFsFromArray(int64_t& index, const double* x_ref);
  void SetNewDelayedAngularDOFsFromArray(int64_t& index, const double* x_ref);

  std::vector<double> GetNewDelayedAngularDOFsAsSTLVector();
  void SetNewDelayedAngularDOFsFromSTLVector(const std::vector<double>& stl_vector);

  std::vector<double> GetOldDelayedAngularDOFsAsSTLVector();
  void SetOldDelayedAngularDOFsFromSTLVector(const std::vector<double>& stl_vector);

  void SetDelayedPsiOld2New();
  void SetDelayedPsiNew2Old();
};

#endif // CHI_ANGLEAGGREGATION_H
