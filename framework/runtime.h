// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "mpicpp-lite/mpicpp-lite.h"
#include "caliper/cali-manager.h"
#include <utility>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory>
#include <filesystem>

namespace mpi = mpicpp_lite;

namespace opensn
{

static const std::string program = "OpenSn";

class FieldFunctionInterpolation;
class Solver;
class MultiGroupXS;
class FieldFunction;
class SpatialDiscretization;
class Timer;
class Logger;
class PostProcessor;
class Object;
class Function;

extern mpi::Communicator mpi_comm;
extern Logger& log;
extern Timer program_timer;
extern bool use_caliper;
extern std::string cali_config;
extern cali::ConfigManager cali_mgr;
extern bool suppress_color;
extern std::filesystem::path input_path;

/// Global stack of handlers
extern std::vector<std::shared_ptr<FieldFunctionInterpolation>> field_func_interpolation_stack;
extern std::vector<std::shared_ptr<MultiGroupXS>> multigroup_xs_stack;
extern std::vector<std::shared_ptr<FieldFunction>> field_function_stack;
extern std::vector<std::shared_ptr<Object>> object_stack;
extern std::vector<std::shared_ptr<SpatialDiscretization>> sdm_stack;
extern std::vector<std::shared_ptr<PostProcessor>> postprocessor_stack;
extern std::vector<std::shared_ptr<Function>> function_stack;

/**
 * Attempts to retrieve an object of base-type `shared_ptr<T>` at the given handle. It then attempts
 * to cast it to type `shared_ptr<R>` and, if  successful, will return a reference of type R&.
 *
 * Example usage:
 *
 * \code
 * const auto& surf_mesh = GetStackItem<SurfaceMesh>(object_stack, surface_hndl);
 * \endcode
 */
template <class R, class T>
static R&
GetStackItem(std::vector<std::shared_ptr<T>>& stack,
             const size_t handle,
             const std::string& calling_function_name = "Unknown")
{
  try
  {
    std::shared_ptr<T>& item = stack.at(handle);
    std::shared_ptr<R> ret_item = std::dynamic_pointer_cast<R>(item);
    if (not ret_item)
      throw std::logic_error("GetStackItem: Invalid return type used. Calling function: " +
                             calling_function_name);
    return *ret_item;
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range("GetStackItem: Invalid handle used. Calling function: " +
                            calling_function_name);
  }
}

/**
 * Attempts to obtain object of type `shared_ptr<T>` at the given
 * handle of a stack with parent type P.
 *
 * Example usage:
 *
 * \code
 * auto surf_mesh_ptr =
 *   GetStackItemPtrAsType<SurfaceMesh>(object_stack, surf_mesh_hndle, fname);
 * \endcode
 */
template <class T, class P>
static std::shared_ptr<T>
GetStackItemPtrAsType(std::vector<std::shared_ptr<P>>& stack,
                      const size_t handle,
                      const std::string& calling_function_name = "Unknown")
{
  std::shared_ptr<P> item_type_P;
  try
  {
    item_type_P = stack.at(handle);
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range("GetStackItem: Invalid handle used. Calling function: " +
                            calling_function_name);
  }

  auto item_type_T = std::dynamic_pointer_cast<T>(item_type_P);
  if (not item_type_T)
    throw std::logic_error(calling_function_name + " Failed to cast to requested type");

  return item_type_T;
}

/**
 * Attempts to obtain object of type `shared_ptr<T>` at the given
 * handle of a stack ALSO OF TYPE T.
 *
 * Example usage:
 *
 * \code
 * auto surf_mesh_ptr = GetStackItemPtr(object_stack, surf_mesh_hndle, fname);
 * \endcode
 */
template <class T>
static std::shared_ptr<T>&
GetStackItemPtr(std::vector<std::shared_ptr<T>>& stack,
                const size_t handle,
                const std::string& calling_function_name = "Unknown")
{
  try
  {
    std::shared_ptr<T>& item = stack.at(handle);
    return item;
  }
  catch (const std::out_of_range& oor)
  {
    throw std::out_of_range("GetStackItem: Invalid handle used. Calling function: " +
                            calling_function_name);
  }
}

/// Initializes all necessary items
int Initialize();

/// Finalize the run
void Finalize();

/// Gets the version string.
std::string GetVersionStr();

} // namespace opensn
