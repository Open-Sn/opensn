#include "sweepscheduler.h"

using namespace chi_mesh::sweep_management;

//=============================================== Phi functions
/**Sets the location where flux moments are to be written.*/
void SweepScheduler::SetDestinationPhi(std::vector<double> &in_destination_phi)
{
  sweep_chunk_.SetDestinationPhi(in_destination_phi);
}

/**Sets all elements of the output vector to zero.*/
void SweepScheduler::ZeroDestinationPhi()
{ sweep_chunk_.ZeroDestinationPhi();
}

/**Returns a reference to the output flux moments vector.*/
std::vector<double>& SweepScheduler::GetDestinationPhi()
{
  return sweep_chunk_.GetDestinationPhi();
}


//=============================================== Psi functions
/**Sets the location where angular fluxes are to be written.*/
void SweepScheduler::SetDestinationPsi(std::vector<double>& in_destination_psi)
{
  sweep_chunk_.SetDestinationPsi(in_destination_psi);
}

/**Sets all elements of the output angular flux vector to zero.*/
void SweepScheduler::ZeroDestinationPsi()
{ sweep_chunk_.ZeroDestinationPsi();
}

/**Returns a reference to the output angular flux vector.*/
std::vector<double>& SweepScheduler::GetDestinationPsi()
{
  return sweep_chunk_.GetDestinationPsi();
}


/** Resets all the incoming intra-location and inter-location
   * cyclic interfaces.*/
void SweepScheduler::ZeroIncomingDelayedPsi()
{
  angle_agg_.ZeroIncomingDelayedPsi();
}

/** Resets all the outgoing intra-location and inter-location
 * cyclic interfaces.*/
void SweepScheduler::ZeroOutgoingDelayedPsi()
{
  angle_agg_.ZeroOutgoingDelayedPsi();
}

/** Clear the output angular flux vector, the flux moments
   * vector, and the outgoing delayed psi.
   */
void SweepScheduler::ZeroOutputFluxDataStructures()
{
  ZeroDestinationPsi();
  ZeroDestinationPhi();
  ZeroOutgoingDelayedPsi();
}

/**Activates or deactives the surface src flag.*/
void SweepScheduler::SetBoundarySourceActiveFlag(bool flag_value)
{
  sweep_chunk_.SetBoundarySourceActiveFlag(flag_value);
}
