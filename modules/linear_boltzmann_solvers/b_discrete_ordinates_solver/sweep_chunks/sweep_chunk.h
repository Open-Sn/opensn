#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/lbs_structs.h"
#include "modules/linear_boltzmann_solvers/a_lbs_solver/groupset/lbs_groupset.h"
#include <functional>

namespace opensn
{
namespace lbs
{

/**Sweep work function*/
class SweepChunk
{
public:
  /**
   * Convenient typdef for the moment call back function. See moment_callbacks.
   *  Arguments are:
   *  SweepChunk *
   *  AngleSet *
   */
  typedef std::function<void(SweepChunk* sweeper, AngleSet* angle_set)> MomentCallbackF;

  /**
   * Functions of type MomentCallbackF can be added to the moment_callbacks
   * vector and these can be called from within functions taking a
   * LBSGroupset instance. The intention is that this function can
   * be used as a general interface to retrieve angular flux values
   */
  std::vector<MomentCallbackF> moment_callbacks;

  SweepChunk(std::vector<double>& in_destination_phi,
             std::vector<double>& in_destination_psi,
             const MeshContinuum& grid,
             const SpatialDiscretization& discretization,
             const std::vector<lbs::UnitCellMatrices>& unit_cell_matrices,
             std::vector<lbs::CellLBSView>& cell_transport_views,
             const std::vector<double>& source_moments,
             const lbs::LBSGroupset& groupset,
             const std::map<int, std::shared_ptr<MultiGroupXS>>& xs,
             int num_moments,
             int max_num_cell_dofs)
    : grid_(grid),
      discretization_(discretization),
      unit_cell_matrices_(unit_cell_matrices),
      cell_transport_views_(cell_transport_views),
      source_moments_(source_moments),
      groupset_(groupset),
      xs_(xs),
      num_moments_(num_moments),
      max_num_cell_dofs_(max_num_cell_dofs),
      save_angular_flux_(not in_destination_psi.empty()),
      groupset_angle_group_stride_(groupset_.psi_uk_man_.NumberOfUnknowns() *
                                   groupset_.groups_.size()),
      groupset_group_stride_(groupset_.groups_.size()),
      destination_phi(&in_destination_phi),
      destination_psi(&in_destination_psi)
  {
  }

  /**Sweep chunks should override this.*/
  virtual void Sweep(AngleSet& angle_set) {}

  /**Sets the currently active FLUx Data Structure*/
  virtual void SetAngleSet(AngleSet& angle_set) {}

  /**For cell-by-cell methods or computing the residual on a single cell.*/
  virtual void SetCell(Cell const* cell_ptr, AngleSet& angle_set) {}

  virtual ~SweepChunk() = default;

protected:
  friend class SweepScheduler;
  /**Sets the location where flux moments are to be written.*/
  void SetDestinationPhi(std::vector<double>& in_destination_phi)
  {
    destination_phi = (&in_destination_phi);
  }

  /**Sets all elements of the output vector to zero.*/
  void ZeroDestinationPhi() { (*destination_phi).assign((*destination_phi).size(), 0.0); }

  /**Returns a reference to the output flux moments vector.*/
  std::vector<double>& GetDestinationPhi() { return *destination_phi; }

  /**Sets the location where angular fluxes are to be written.*/
  void SetDestinationPsi(std::vector<double>& in_destination_psi)
  {
    destination_psi = (&in_destination_psi);
  }

  /**Sets all elements of the output angular flux vector to zero.*/
  void ZeroDestinationPsi() { (*destination_psi).assign((*destination_psi).size(), 0.0); }

  /**Returns a reference to the output angular flux vector.*/
  std::vector<double>& GetDestinationPsi() { return *destination_psi; }

  /**Activates or deactives the surface src flag.*/
  void SetBoundarySourceActiveFlag(bool flag_value) // Done
  {
    surface_source_active = flag_value;
  }

  /**Returns the surface src-active flag.*/
  bool IsSurfaceSourceActive() const { return surface_source_active; }

  const MeshContinuum& grid_;
  const SpatialDiscretization& discretization_;
  const std::vector<lbs::UnitCellMatrices>& unit_cell_matrices_;
  std::vector<lbs::CellLBSView>& cell_transport_views_;
  const std::vector<double>& source_moments_;
  const lbs::LBSGroupset& groupset_;
  const std::map<int, std::shared_ptr<MultiGroupXS>>& xs_;
  const int num_moments_;
  const int max_num_cell_dofs_;
  const bool save_angular_flux_;
  const size_t groupset_angle_group_stride_;
  const size_t groupset_group_stride_;

private:
  std::vector<double>* destination_phi;
  std::vector<double>* destination_psi;
  bool surface_source_active = false;
};

} // namespace lbs
} // namespace opensn
