#include "framework/object.h"

namespace opensn
{

InputParameters
Object::GetInputParameters()
{
  return {}; // Returns an empty block
}

Object::Object()
{
}

Object::Object(const InputParameters&)
{
}

void
Object::SetStackID(size_t stack_id)
{
  stack_id_ = stack_id;
}

size_t
Object::StackID() const
{
  return stack_id_;
}

void
Object::PushOntoStack(std::shared_ptr<Object>& new_object)
{
  object_stack.push_back(new_object);
  new_object->SetStackID(object_stack.size() - 1);
}

} // namespace opensn
