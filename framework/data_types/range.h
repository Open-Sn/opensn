// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>
#include <vector>

namespace opensn
{
/**
 * Returns a range of number according to the logic of the parameters.
 *
 * \param start First number in the sequence.
 * \param end   Termination criteria. If the delta is positive then
 *              the sequence will terminate if i>=end, otherwise if the
 *              delta is negative the sequence will terminate if i<=end
 * \param delta Cannot be 0. Default 1. Can be negative.*/
template <typename T, typename D = int>
std::vector<T>
Range(T start, T end, D delta = 1)
{
  static_assert(std::is_signed<D>::value, "Range delta parameter must be signed");
  const bool forward = (delta > 0);

  std::vector<T> sequence = {};

  if (forward and start >= end)
    return sequence;
  if (not forward and start <= end)
    return sequence;

  T i = start;
  bool terminate = false;
  while (not terminate)
  {
    sequence.push_back(i);
    i += delta;

    if (forward and i >= end)
      terminate = true;
    if (not forward and i <= end)
      terminate = true;

    // Wrap-around check
    if (not forward and i > start)
      terminate = true;
  }

  return sequence;
}

} // namespace opensn
