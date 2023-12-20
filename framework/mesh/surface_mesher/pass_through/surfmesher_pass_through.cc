#include "framework/mesh/surface_mesher/pass_through/surfmesher_pass_through.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

namespace opensn
{

void
SurfaceMesherPassthrough::Execute()
{
  Chi::log.Log0Verbose1() << "SurfaceMesherPassthrough executed";
}

} // namespace opensn
