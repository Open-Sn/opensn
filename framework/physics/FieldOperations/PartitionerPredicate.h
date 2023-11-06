#pragma once

#include "framework/physics/FieldOperations/field_operation.h"

namespace chi
{
class GraphPartitioner;
}

namespace chi_physics::field_operations
{

class PartitionerPredicate : public FieldOperation
{
public:
  static chi::InputParameters GetInputParameters();
  explicit PartitionerPredicate(const chi::InputParameters& params);

  void Execute() override;

protected:
  chi::GraphPartitioner& partitioner_;
  const chi::ParameterBlock result_field_param_;
  const size_t num_partitions_;
  const size_t result_component_;
};

} // namespace chi_physics::field_operations
