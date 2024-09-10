// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

extern "C"
{
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"
}
#include "LuaBridge/LuaBridge.h"
#include "framework/logging/log.h"
#include "framework/field_functions/interpolation/ffinterpolation.h"
#include "framework/math/quadratures/angular/product_quadrature.h"

namespace luabridge
{

template <>
struct Stack<opensn::Logger::LOG_LVL> : Enum<opensn::Logger::LOG_LVL,
                                             opensn::Logger::LOG_0,
                                             opensn::Logger::LOG_0WARNING,
                                             opensn::Logger::LOG_0ERROR,
                                             opensn::Logger::LOG_0VERBOSE_0,
                                             opensn::Logger::LOG_0VERBOSE_1,
                                             opensn::Logger::LOG_0VERBOSE_2,
                                             opensn::Logger::LOG_ALL,
                                             opensn::Logger::LOG_ALLWARNING,
                                             opensn::Logger::LOG_ALLERROR,
                                             opensn::Logger::LOG_ALLVERBOSE_0,
                                             opensn::Logger::LOG_ALLVERBOSE_1,
                                             opensn::Logger::LOG_ALLVERBOSE_2>
{
};

template <>
struct Stack<opensn::FieldFunctionInterpolationOperation>
  : Enum<opensn::FieldFunctionInterpolationOperation,
         opensn::FieldFunctionInterpolationOperation::OP_SUM,
         opensn::FieldFunctionInterpolationOperation::OP_AVG,
         opensn::FieldFunctionInterpolationOperation::OP_MAX,
         opensn::FieldFunctionInterpolationOperation::OP_SUM_FUNC,
         opensn::FieldFunctionInterpolationOperation::OP_AVG_FUNC,
         opensn::FieldFunctionInterpolationOperation::OP_MAX_FUNC>
{
};

template <>
struct Stack<opensn::ProductQuadratureType>
  : Enum<opensn::ProductQuadratureType,
         opensn::ProductQuadratureType::UNKNOWN,
         opensn::ProductQuadratureType::GAUSS_LEGENDRE,
         opensn::ProductQuadratureType::GAUSS_CHEBYSHEV,
         opensn::ProductQuadratureType::GAUSS_LEGENDRE_CHEBYSHEV,
         opensn::ProductQuadratureType::CUSTOM_QUADRATURE>
{
};

} // namespace luabridge
