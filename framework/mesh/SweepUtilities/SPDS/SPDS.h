#pragma once

#include "framework/mesh/SweepUtilities/SPLS/SPLS.h"
#include "framework/mesh/chi_mesh.h"

#include <memory>

namespace chi_mesh::sweep_management
{

/**Contains multiple levels*/
class SPDS
{
public:
  SPDS(const chi_mesh::Vector3& in_omega, const chi_mesh::MeshContinuum& in_grid, bool verbose)
    : omega_(in_omega), grid_(in_grid), verbose_(verbose)
  {
  }

  const chi_mesh::MeshContinuum& Grid() const { return grid_; }
  const chi_mesh::Vector3& Omega() const { return omega_; }
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
  chi_mesh::Vector3 omega_;

  const chi_mesh::MeshContinuum& grid_;

  SPLS spls_;

  std::vector<int> location_dependencies_;
  std::vector<int> location_successors_;
  std::vector<int> delayed_location_dependencies_;
  std::vector<int> delayed_location_successors_;

  std::vector<std::pair<int, int>> local_cyclic_dependencies_;

  std::vector<std::vector<FaceOrientation>> cell_face_orientations_;

  bool verbose_ = false;

  /**Populates cell relationships and cell_face_orientations.*/
  void PopulateCellRelationships(const chi_mesh::Vector3& omega,
                                 std::set<int>& location_dependencies,
                                 std::set<int>& location_successors,
                                 std::vector<std::set<std::pair<int, double>>>& cell_successors);

  void PrintedGhostedGraph() const;
};

} // namespace chi_mesh::sweep_management
