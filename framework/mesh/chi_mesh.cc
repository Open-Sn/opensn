#include "chi_mesh.h"
#include "MeshHandler/chi_meshhandler.h"
#include "chi_runtime.h"

namespace chi_mesh
{

MeshHandler&
GetCurrentHandler()
{
  if (Chi::meshhandler_stack.empty())
    throw std::logic_error("chi_mesh::GetCurrentHandler: No handlers on stack");

  return Chi::GetStackItem<MeshHandler>(Chi::meshhandler_stack, Chi::current_mesh_handler);
}

size_t
PushNewHandlerAndGetIndex()
{
  Chi::meshhandler_stack.push_back(std::make_shared<MeshHandler>());

  int index = (int)Chi::meshhandler_stack.size() - 1;
  Chi::current_mesh_handler = index;

  return index;
}

} // namespace chi_mesh
