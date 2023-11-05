#pragma once

#include <petscksp.h>

#include "opensn/modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_context.h"

namespace lbs
{
int WGDSA_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output);
int WGDSA_TGDSA_PreConditionerMult2(lbs::WGSContext<Mat, Vec, KSP>& gs_context_ptr,
                                    Vec phi_input,
                                    Vec pc_output);
int MIP_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output);
} // namespace lbs
