#include "chi_mesh.h"
#include<iostream>

#include "MeshHandler/chi_meshhandler.h"
#include "chi_runtime.h"


//###################################################################
/**Obtains a reference to the current mesh handler from the global stack.
 *
 * If the stack is empty this routine will through `std::logic_error`.
\author Jan*/
chi_mesh::MeshHandler& chi_mesh::GetCurrentHandler()
{
  if (Chi::meshhandler_stack.empty())
    throw std::logic_error("chi_mesh::GetCurrentHandler: No handlers on stack");

  return Chi::GetStackItem<chi_mesh::MeshHandler>(Chi::meshhandler_stack,
                                                  Chi::current_mesh_handler);
}

//###################################################################
/**Adds a new mesh handler to the stack, sets it as the current handler
 * and returns a handle to it.*/
size_t chi_mesh::PushNewHandlerAndGetIndex()
{
  Chi::meshhandler_stack.push_back(std::make_shared<chi_mesh::MeshHandler>());

  int index = (int)Chi::meshhandler_stack.size() - 1;
  Chi::current_mesh_handler = index;

  return index;
}


