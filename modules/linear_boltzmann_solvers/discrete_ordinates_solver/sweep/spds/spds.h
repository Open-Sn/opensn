#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/sweep.h"
#include <memory>

namespace opensn
{
namespace lbs
{

/**Contains multiple levels*/
class SPDS
{
public:
  SPDS(const Vector3& omega, const MeshContinuum& grid, bool verbose)
    : omega_(omega), grid_(grid), verbose_(verbose)
  {
  }

  const MeshContinuum& Grid() const { return grid_; }
  const Vector3& Omega() const { return omega_; }
  const SPLS& GetSPLS() const { return spls_; }

  typedef std::vector<int> VecInt;
  const VecInt& GetLocationDependencies() const { return location_dependencies_; }
  const VecInt& GetLocationSuccessors() const { return location_successors_; }
  const VecInt& GetDelayedLocationDependencies() const { return delayed_location_dependencies_; }
  const VecInt& GetDelayedLocationSuccessors() const { return delayed_location_successors_; }
  const std::vector<std::pair<int, int>>& GetLocalCyclicDependencies() const
  {
    return local_cyclic_dependencies_;
  }
  const std::vector<std::vector<FaceOrientation>>& CellFaceOrientations() const
  {
    return cell_face_orientations_;
  }

  /** Given a location J index, maps to a predecessor location.*/
  int MapLocJToPrelocI(int locJ) const;
  /** Given a location J index, maps to a dependent location.*/
  int MapLocJToDeplocI(int locJ) const;

  virtual ~SPDS() = default;

protected:
  Vector3 omega_;

  const MeshContinuum& grid_;

  SPLS spls_;

  std::vector<int> location_dependencies_;
  std::vector<int> location_successors_;
  std::vector<int> delayed_location_dependencies_;
  std::vector<int> delayed_location_successors_;

  std::vector<std::pair<int, int>> local_cyclic_dependencies_;

  std::vector<std::vector<FaceOrientation>> cell_face_orientations_;

  bool verbose_ = false;

  /**Populates cell relationships and cell_face_orientations.*/
  void PopulateCellRelationships(const Vector3& omega,
                                 std::set<int>& location_dependencies,
                                 std::set<int>& location_successors,
                                 std::vector<std::set<std::pair<int, double>>>& cell_successors);

  void PrintedGhostedGraph() const;
};

} // namespace lbs
} // namespace opensn
