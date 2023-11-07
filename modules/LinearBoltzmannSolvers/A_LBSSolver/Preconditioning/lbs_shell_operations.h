#pragma once

#include <petscksp.h>

#include "modules/LinearBoltzmannSolvers/A_LBSSolver/IterativeMethods/wgs_context.h"

namespace lbs
{
/**Applies WGDSA or TGDSA to the given input vector.*/
int WGDSA_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output);
/**Applies WGDSA or TGDSA to the given input vector.*/
int WGDSA_TGDSA_PreConditionerMult2(lbs::WGSContext<Mat, Vec, KSP>& gs_context_ptr,
                                    Vec phi_input,
                                    Vec pc_output);
/**Applies TGDSA to the given input vector.*/
int MIP_TGDSA_PreConditionerMult(PC pc, Vec phi_input, Vec pc_output);
} // namespace lbs
