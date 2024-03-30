#pragma once

#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep/angle_aggregation/angle_aggregation.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_solver/sweep_chunks/sweep_chunk.h"

namespace opensn
{
namespace lbs
{

class SweepChunk;

enum class SchedulingAlgorithm
{
  FIRST_IN_FIRST_OUT = 1, ///< FIFO
  DEPTH_OF_GRAPH = 2      ///< DOG
};

typedef AngleSetGroup TAngleSetGroup;
typedef AngleSet TAngleSet;
typedef STDG TGSPO;

typedef std::vector<TGSPO> TLEVELED_GRAPH;

class SweepScheduler
{
private:
  SchedulingAlgorithm scheduler_type_;
  AngleAggregation& angle_agg_;

  struct RULE_VALUES
  {
    std::shared_ptr<TAngleSet> angle_set;
    int depth_of_graph;
    int sign_of_omegax;
    int sign_of_omegay;
    int sign_of_omegaz;
    size_t set_index;

    explicit RULE_VALUES(std::shared_ptr<TAngleSet>& ref_as) : angle_set(ref_as)
    {
      depth_of_graph = 0;
      set_index = 0;
      sign_of_omegax = 1;
      sign_of_omegay = 1;
      sign_of_omegaz = 1;
    }
  };
  std::vector<RULE_VALUES> rule_values_;

  SweepChunk& sweep_chunk_;

public:
  SweepScheduler(SchedulingAlgorithm scheduler_type,
                 AngleAggregation& angle_agg,
                 SweepChunk& sweep_chunk);

  AngleAggregation& AngleAgg() { return angle_agg_; }

  /**
   * This is the entry point for sweeping.
   */
  void Sweep();

  /**
   * Returns the referenced sweep chunk.
   */
  SweepChunk& GetSweepChunk();

private:
  /**
   * Applies a First-In-First-Out sweep scheduling.
   */
  void ScheduleAlgoFIFO(SweepChunk& sweep_chunk);

  /**
   * Initializes the Depth-Of-Graph algorithm.
   */
  void InitializeAlgoDOG();

  /**
   * Executes the Depth-Of-Graph algorithm.
   */
  void ScheduleAlgoDOG(SweepChunk& sweep_chunk);

public:
  /**
   * Sets the location where flux moments are to be written.
   */
  void SetDestinationPhi(std::vector<double>& destination_phi);

  /**
   * Sets all elements of the output vector to zero.
   */
  void ZeroDestinationPhi();

  /**
   * Returns a reference to the output flux moments vector.
   */
  std::vector<double>& GetDestinationPhi();

  /**
   * Sets the location where angular fluxes are to be written.
   */
  void SetDestinationPsi(std::vector<double>& destination_psi);

  /**
   * Sets all elements of the output angular flux vector to zero.
   */
  void ZeroDestinationPsi();

  /**
   * Returns a reference to the output angular flux vector.
   */
  std::vector<double>& GetDestinationPsi();

  /** Resets all the incoming intra-location and inter-location
   * cyclic interfaces.
   */
  void ZeroIncomingDelayedPsi();

  /** Resets all the outgoing intra-location and inter-location
   * cyclic interfaces.
   */
  void ZeroOutgoingDelayedPsi();

  /**
   * Clear the output angular flux vector, the flux moments vector, and the outgoing delayed psi.
   */
  void ZeroOutputFluxDataStructures();

  /**
   * Activates or deactives the surface src flag.
   */
  void SetBoundarySourceActiveFlag(bool flag_value);
};

} // namespace lbs
} // namespace opensn
