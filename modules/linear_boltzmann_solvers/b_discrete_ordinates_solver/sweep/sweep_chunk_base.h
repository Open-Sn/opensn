#pragma once

#include "modules/linear_boltzmann_solvers/b_discrete_ordinates_solver/sweep/angle_aggregation/angle_aggregation.h"

#include <functional>

namespace opensn
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

  SweepChunk(std::vector<double>& in_destination_phi, std::vector<double>& in_destination_psi)
    : destination_phi(&in_destination_phi), destination_psi(&in_destination_psi)
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

protected:
  /**Returns the surface src-active flag.*/
  bool IsSurfaceSourceActive() const { return surface_source_active; }

private:
  std::vector<double>* destination_phi;
  std::vector<double>* destination_psi;
  bool surface_source_active = false;
};

} // namespace opensn
