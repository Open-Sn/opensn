#pragma once

#include "opensn/framework/chi_runtime.h"
#include "opensn/framework/parameters/input_parameters.h"

class ChiObject
{
private:
  size_t stack_id_ = Chi::SIZE_T_INVALID;

public:
  /**Returns the input parameters. For the base ChiObject, there
   * are now parameters loaded.*/
  static chi::InputParameters GetInputParameters();

  /**Default constructor. This will be removed in future.*/
  ChiObject();

  /**Constructor with input parameters.*/
  explicit ChiObject(const chi::InputParameters& params);

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
  virtual void PushOntoStack(std::shared_ptr<ChiObject>& new_object);

  virtual ~ChiObject() = default;
};
