#pragma once

#include "framework/runtime.h"
#include "framework/parameters/input_parameters.h"

namespace opensn
{

class Object
{
private:
  size_t stack_id_ = Chi::SIZE_T_INVALID;

public:
  /**Returns the input parameters. For the base Object, there
   * are now parameters loaded.*/
  static InputParameters GetInputParameters();

  /**Default constructor. This will be removed in future.*/
  Object();

  /**Constructor with input parameters.*/
  explicit Object(const InputParameters& params);

  // Setters
  /**Sets the stack id of the object. This allows this
   * object to know its place in the global space.*/
  void SetStackID(size_t stack_id);

  // Getters
  /**Returns the stack id of this object. This can be used
   * with input language to connect objects together.*/
  size_t StackID() const;

  /**An overridable callback that is called by the ObjectMaker and by default
   * adds the object onto the object stack. This function can be used to
   * place the object on a different stack.*/
  virtual void PushOntoStack(std::shared_ptr<Object>& new_object);

  virtual ~Object() = default;
};

} // namespace opensn
