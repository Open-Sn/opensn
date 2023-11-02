#include "lbs_solver.h"

#include "chi_runtime.h"
#include "chi_log.h"

//###################################################################
/** Initialize the solver.*/
void lbs::LBSSolver::Initialize()
{
  PerformInputChecks();                //a assigns num_groups and grid
  PrintSimHeader();                    //b

  MPI_Barrier(Chi::mpi.comm);

  InitMaterials();                     //c
  InitializeSpatialDiscretization();   //d
  InitializeGroupsets();               //e
  ComputeNumberOfMoments();            //f
  InitializeParrays();                 //g
  InitializeBoundaries();              //h
  InitializePointSources();            //i

  source_event_tag_ = Chi::log.GetRepeatingEventTag("Set Source");
}
