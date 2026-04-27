// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <iomanip>
#include <ostream>
#include <string>

namespace opensn
{

enum class NumericFormatKind
{
  FIXED,
  SCIENTIFIC,
};

struct NumericFormat
{
  int precision = 6;
  NumericFormatKind kind = NumericFormatKind::SCIENTIFIC;
};

inline NumericFormat
Fixed(const int precision)
{
  return {precision, NumericFormatKind::FIXED};
}

inline NumericFormat
Scientific(const int precision)
{
  return {precision, NumericFormatKind::SCIENTIFIC};
}

inline void
AppendNumericField(std::ostream& out,
                   const std::string& label,
                   const std::size_t value,
                   const bool leading_separator = true)
{
  out << (leading_separator ? ", " : " ") << label << " = " << value;
}

inline void
AppendNumericField(std::ostream& out,
                   const std::string& label,
                   const double value,
                   const NumericFormat format,
                   const bool leading_separator = true)
{
  out << (leading_separator ? ", " : " ") << label << " = ";
  if (format.kind == NumericFormatKind::FIXED)
    out << std::fixed;
  else
    out << std::scientific;
  out << std::setprecision(format.precision) << value << std::defaultfloat;
}

} // namespace opensn
