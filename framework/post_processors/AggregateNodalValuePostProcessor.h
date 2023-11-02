#ifndef CHITECH_AGGREGATENODALVALUEPOSTPROCESSOR_H
#define CHITECH_AGGREGATENODALVALUEPOSTPROCESSOR_H

#include "PostProcessor.h"
#include "physics/FieldFunction/GridBasedFieldFunctionInterface.h"
#include "mesh/LogicalVolume/LogicalVolumeInterface.h"

namespace chi_mesh
{
class LogicalVolume;
}
namespace chi_physics
{
class FieldFunctionGridBased;
}

namespace chi
{

class AggregateNodalValuePostProcessor
  : public PostProcessor,
    public chi_physics::GridBasedFieldFunctionInterface,
    public chi_mesh::LogicalVolumeInterface
{
public:
  static InputParameters GetInputParameters();
  explicit AggregateNodalValuePostProcessor(const InputParameters& params);

  void Execute(const Event& event_context) override;

protected:
  void Initialize();

  const std::string operation_;
  bool initialized_ = false;
  std::vector<uint64_t> cell_local_ids_;

};

} // namespace chi

#endif // CHITECH_AGGREGATENODALVALUEPOSTPROCESSOR_H
