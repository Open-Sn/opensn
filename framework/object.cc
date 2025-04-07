// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/object.h"
#include <limits>

namespace opensn
{

InputParameters
Object::GetInputParameters()
{
  return {}; // Returns an empty block
}

Object::Object() : stack_id_(std::numeric_limits<std::size_t>::max())
{
}

Object::Object(const InputParameters&) : stack_id_(std::numeric_limits<std::size_t>::max())
{
}

void
Object::SetStackID(size_t stack_id)
{
  stack_id_ = stack_id;
}

size_t
Object::GetStackID() const
{
  return stack_id_;
}

void
Object::PushOntoStack(std::shared_ptr<Object> new_object)
{
  object_stack.push_back(new_object);
  new_object->SetStackID(object_stack.size() - 1);
}

} // namespace opensn
