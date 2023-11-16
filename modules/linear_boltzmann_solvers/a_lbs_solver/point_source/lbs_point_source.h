#pragma once

#include <utility>

#include "framework/mesh/mesh.h"

namespace opensn
{
namespace lbs
{

class PointSource
{
private:
  const Vector3 location_;
  const std::vector<double> groupwise_strength_;

public:
  struct ContainingCellInfo
  {
    double volume_weight;
    uint64_t cell_local_id;
    std::vector<double> shape_values;
    std::vector<double> node_weights;
  };

private:
  std::vector<ContainingCellInfo> m_containing_cells_;

public:
  PointSource(const Vector3& location, std::vector<double> strength)
    : location_(location), groupwise_strength_(std::move(strength))
  {
  }

  const Vector3& Location() const { return location_; }

  const std::vector<double>& Strength() const { return groupwise_strength_; }

  void AddContainingCellInfo(double volume_weight,
                             uint64_t cell_local_id,
                             std::vector<double> shape_values,
                             std::vector<double> node_weights)
  {
    m_containing_cells_.push_back(ContainingCellInfo{
      volume_weight, cell_local_id, std::move(shape_values), std::move(node_weights)});
  }

  const std::vector<ContainingCellInfo>& ContainingCellsInfo() const { return m_containing_cells_; }

  void ClearInitializedInfo() { m_containing_cells_.clear(); }
};

} // namespace lbs
} // namespace opensn
