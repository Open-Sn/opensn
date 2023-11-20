#pragma once

#include <utility>
#include <vector>
#include <string>
#include <stdexcept>
#include <memory>

#include "framework/mpi/mpi.h"

namespace opensn
{

class MeshHandler;
typedef std::shared_ptr<MeshHandler> MeshHandlerPtr;

class SurfaceMesh;
typedef std::shared_ptr<SurfaceMesh> SurfaceMeshPtr;

class FieldFunctionInterpolation;
typedef FieldFunctionInterpolation FFInterp;
typedef std::shared_ptr<FFInterp> FFInterpPtr;

class UnpartitionedMesh;
typedef std::shared_ptr<UnpartitionedMesh> UnpartitionedMeshPtr;
typedef UnpartitionedMeshPtr UnpartMeshPtr;

class Solver;
class Material;
class MultiGroupXS;
class FieldFunction;
class Function;

typedef std::shared_ptr<Material> MaterialPtr;
typedef std::shared_ptr<MultiGroupXS> MultiGroupXSPtr;
typedef std::shared_ptr<FieldFunction> FieldFunctionPtr;
typedef std::shared_ptr<Function> FunctionPtr;

class AngularQuadrature;
class SpatialDiscretization;

typedef std::shared_ptr<AngularQuadrature> AngularQuadraturePtr;
typedef std::shared_ptr<SpatialDiscretization> SpatialDiscretizationPtr;

class UnknownManager;

class Timer;
class Logger;
class PostProcessor;
typedef std::shared_ptr<PostProcessor> PostProcessorPtr;

class Object;
typedef std::shared_ptr<Object> ChiObjectPtr;

extern MPI_Info& mpi;
extern Logger& log;
extern Timer program_timer;

/**General utilities in ChiTech*/
class Chi
{
public:
  /** Global stack of handlers */
  static std::vector<MeshHandlerPtr> meshhandler_stack;
  static int current_mesh_handler;

  static std::vector<SurfaceMeshPtr> surface_mesh_stack;
  static std::vector<FFInterpPtr> field_func_interpolation_stack;
  static std::vector<UnpartMeshPtr> unpartitionedmesh_stack;

  static std::vector<MaterialPtr> material_stack;
  static std::vector<MultiGroupXSPtr> multigroup_xs_stack;
  static std::vector<FieldFunctionPtr> field_function_stack;

  static std::vector<AngularQuadraturePtr> angular_quadrature_stack;

  static std::vector<ChiObjectPtr> object_stack;
  static std::vector<SpatialDiscretizationPtr> sdm_stack;
  static std::vector<PostProcessorPtr> postprocessor_stack;
  static std::vector<FunctionPtr> function_stack;

  static const size_t SIZE_T_INVALID = ((size_t)-1);

  static bool suppress_color_;

  /**Customized exceptions.*/
  class RecoverableException : public std::runtime_error
  {
  public:
    explicit RecoverableException(const char* message)
      : std::runtime_error(std::string("RecoverableException: ") + std::string(message))
    {
    }
    explicit RecoverableException(const std::string& message)
      : std::runtime_error(std::string("RecoverableException: ") + message)
    {
    }
    RecoverableException(const std::string& prefix, const std::string& message)
      : std::runtime_error(prefix + message)
    {
    }

    ~RecoverableException() noexcept override = default;
  };

public:
  /// Deleted constructor
  Chi() = delete;
  /// Deleted copy constructor
  Chi(const Chi&) = delete;
  /// Deleted assigment operator
  Chi operator=(const Chi&) = delete;

public:
  /**Attempts to retrieve an object of base-type `shared_ptr<T>` at the given
   * handle. It then attempts to cast it to type `shared_ptr<R>` and, if
   * successful, will return a reference of type R&.
   * \n
   * \n
   * Example usage:
   *
   * \code
   * const auto& surf_mesh = Chi::GetStackItem<SurfaceMesh>(
        Chi::object_stack, surface_hndl);
     // Returns SurfaceMesh&
   * \endcode
   * */
  template <class R, class T>
  static R& GetStackItem(std::vector<std::shared_ptr<T>>& stack,
                         const size_t handle,
                         const std::string& calling_function_name = "Unknown")
  {
    try
    {
      std::shared_ptr<T>& item = stack.at(handle);
      std::shared_ptr<R> ret_item = std::dynamic_pointer_cast<R>(item);
      if (not ret_item)
        throw std::logic_error("chi::GetStackItem: Invalid return type used. "
                               "Calling function: " +
                               calling_function_name);
      return *ret_item;
    }
    catch (const std::out_of_range& oor)
    {
      throw std::out_of_range("chi::GetStackItem: Invalid handle used. "
                              "Calling function: " +
                              calling_function_name);
    }
  }

  /**Attempts to obtain object of type `shared_ptr<T>` at the given
   * handle of a stack with parent type P.
   * \n
   * \n
   * Example usage:
   *
   * \code
   * auto surf_mesh_ptr =
   *   Chi::GetStackItemPtrAsType<SurfaceMesh>(
         Chi::object_stack, surf_mesh_hndle, fname);
     // Returns std::shared_ptr<SurfaceMesh>
   * \endcode
   * */
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
      throw std::out_of_range("chi::GetStackItem: Invalid handle used. "
                              "Calling function: " +
                              calling_function_name);
    }

    auto item_type_T = std::dynamic_pointer_cast<T>(item_type_P);
    if (not item_type_T)
      throw std::logic_error(calling_function_name + "Failed to cast to requested type");

    return item_type_T;
  }

  /**Attempts to obtain object of type `shared_ptr<T>` at the given
   * handle of a stack ALSO OF TYPE T.
   * \n
   * \n
   * Example usage:
   *
   * \code
   * auto surf_mesh_ptr = Chi::GetStackItemPtr(
      Chi::object_stack, surf_mesh_hndle, fname);
     // Returns std::shared_ptr<Object>
   * \endcode
   * */
  template <class T>
  static std::shared_ptr<T>& GetStackItemPtr(std::vector<std::shared_ptr<T>>& stack,
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
      throw std::out_of_range("chi::GetStackItem: Invalid handle used. "
                              "Calling function: " +
                              calling_function_name);
    }
  }
};

/**
 * Initializes all necessary items
 */
int Initialize();

/**
 * Finalize the run
 */
void Finalize();

/**
 * Gets the version string.
 */
std::string GetVersionStr();

/**
 * Exits the program appropriately.
 */
void Exit(int error_code);

} // namespace opensn
