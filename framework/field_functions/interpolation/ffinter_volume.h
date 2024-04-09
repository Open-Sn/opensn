// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/mesh/logical_volume/logical_volume.h"
#include "framework/math/functions/scalar_material_function.h"

#include <petscksp.h>

namespace opensn
{

/**
 * Volume-wise field function interpolation.
 *
 * This interpolator allows the user to obtain quantities by logical
 * volume. If no logical volume is assigned to the method it will
 * default to operating over the entire volume.\n
 * \n
 * The method also supports a few primitive operations:
 *  - OP_VOLUME_AVG. Obtains the volume average of the field function
 *    of interest.
 *  - OP_VOLUME_SUM. Obtains the volume integral of the field function
 *    of interest.
 */
class FieldFunctionInterpolationVolume : public FieldFunctionInterpolation
{
protected:
  std::shared_ptr<LogicalVolume> logical_volume_ = nullptr;
  FieldFunctionInterpolationOperation op_type_ = FieldFunctionInterpolationOperation::OP_SUM;
  double op_value_ = 0.0;
  std::shared_ptr<ScalarMaterialFunction> oper_function_;

private:
  std::vector<uint64_t> cell_local_ids_inside_logvol_;

public:
  FieldFunctionInterpolationVolume()
    : FieldFunctionInterpolation(FieldFunctionInterpolationType::VOLUME)
  {
  }

  virtual ~FieldFunctionInterpolationVolume() {}

  std::shared_ptr<LogicalVolume>& GetLogicalVolume() { return logical_volume_; }

  FieldFunctionInterpolationOperation& GetOperationType() { return op_type_; }

  double& GetOpValue() { return op_value_; }

  void SetOperationFunction(std::shared_ptr<ScalarMaterialFunction> function);

  void Initialize() override;
  void Execute() override;

  std::string GetDefaultFileBaseName() const override { return "ZVFFI"; }
  void ExportPython(std::string base_name) override {}
};

} // namespace opensn
