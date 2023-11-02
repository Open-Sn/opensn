#include "ChiObject.h"

/**Returns the input parameters.*/
chi::InputParameters ChiObject::GetInputParameters()
{
  return {}; // Returns an empty block
}

ChiObject::ChiObject() {}

ChiObject::ChiObject(const chi::InputParameters&) {}

void ChiObject::SetStackID(size_t stack_id) { stack_id_ = stack_id; }

size_t ChiObject::StackID() const { return stack_id_; }

void ChiObject::PushOntoStack(std::shared_ptr<ChiObject>& new_object)
{
  Chi::object_stack.push_back(new_object);
  new_object->SetStackID(Chi::object_stack.size() - 1);
}