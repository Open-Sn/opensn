// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/materials/interpolator/interpolator.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <ranges>

namespace opensn
{

namespace
{

void
ContiguousToNDIndex(std::size_t c_index,
                    std::span<const std::uint32_t> shape,
                    std::span<std::uint32_t> index)
{
  if (index.size() != shape.size())
    throw std::invalid_argument("Index rank does not match shape rank.");

  for (std::size_t d = shape.size(); d-- > 0;)
  {
    if (shape[d] == 0)
      throw std::invalid_argument("Shape entries must be positive.");

    index[d] = static_cast<std::uint32_t>(c_index % shape[d]);
    c_index /= shape[d];
  }

  if (c_index != 0)
    throw std::out_of_range("Flattened index exceeds the shape extent.");
}

bool
PrecursorsAreIdentical(const std::pair<double, std::vector<double>>& lhs,
                       const std::pair<double, std::vector<double>>& rhs)
{
  static auto are_relatively_equal = [](const double lhs_value, const double rhs_value)
  {
    constexpr double tolerance = 1.0e-12;
    if (lhs_value == rhs_value)
      return true;

    double scale = std::min(std::abs(lhs_value), std::abs(rhs_value));
    return scale > 0.0 and std::abs(lhs_value - rhs_value) / scale < tolerance;
  };

  if (not are_relatively_equal(lhs.first, rhs.first) or lhs.second.size() != rhs.second.size())
  {
    return false;
  }

  return std::equal(lhs.second.begin(), lhs.second.end(), rhs.second.begin(), are_relatively_equal);
}

double
GetFractionalYield(const std::vector<MultiGroupXS::Precursor>& precursors,
                   const std::pair<double, std::vector<double>>& precursor_prop)
{
  const auto matches_precursor = [&precursor_prop](const MultiGroupXS::Precursor& precursor)
  {
    return PrecursorsAreIdentical({precursor.decay_constant, precursor.emission_spectrum},
                                  precursor_prop);
  };
  const auto precursor = std::ranges::find_if(precursors, matches_precursor);
  return precursor == precursors.end() ? 0.0 : precursor->fractional_yield;
}

SparseMatrix
Reconstruct(const SparseMatrixProperties& prop, const double*& data)
{
  SparseMatrix result(prop.row_size, prop.col_size);

  auto extract_col_idx = [&prop](std::size_t x) -> std::vector<std::size_t>
  {
    auto view = prop.present_indices |
                std::views::filter([x](const auto& pair) { return pair.first == x; }) |
                std::views::transform([](const auto& pair) { return pair.second; });
    return {std::ranges::begin(view), std::ranges::end(view)};
  };

  const double* end_data;
  for (std::size_t row_idx = 0; row_idx < prop.row_size; ++row_idx)
  {
    result.rowI_indices[row_idx] = extract_col_idx(row_idx);
    end_data = data + result.rowI_indices[row_idx].size();
    result.rowI_values[row_idx] = std::vector<double>(data, end_data);
    data = end_data;
  }
  return result;
}

} // namespace

SparseMatrixProperties&
SparseMatrixProperties::operator+=(const SparseMatrixProperties& prop)
{
  if ((row_size != prop.row_size) or (col_size != prop.col_size))
    throw std::runtime_error("Shape mismatch.");
  present_indices.insert(prop.present_indices.begin(), prop.present_indices.end());
  return *this;
}

SparseMatrixProperties::SparseMatrixProperties(const SparseMatrix& mat)
  : row_size(mat.GetNumRows()), col_size(mat.GetNumCols())
{
  for (std::size_t row_idx = 0; row_idx < mat.GetNumRows(); ++row_idx)
  {
    for (const auto& col_idx : mat.rowI_indices[row_idx])
      present_indices.insert({row_idx, col_idx});
  }
}

XSProperties::XSProperties(const MultiGroupXS& xs)
  : num_groups(xs.GetNumGroups()),
    scattering_order(xs.GetScatteringOrder()),
    num_precursors(xs.GetNumPrecursors()),
    is_fissionable(xs.IsFissionable()),
    adjoint_mode(xs.GetAdjointMode()),
    energy_deposition(xs.GetEnergyDeposition())
{
  for (const auto& transfer_matrix : xs.GetTransferMatrices())
    transfer_matrix_props.emplace_back(transfer_matrix);
  transfer_matrix_props.shrink_to_fit();

  for (const auto& precursor : xs.GetPrecursors())
    InsertPrecursor({precursor.decay_constant, precursor.emission_spectrum});
}

bool
XSProperties::operator==(const XSProperties& other) const
{
  return num_groups == other.num_groups and scattering_order == other.scattering_order and
         num_precursors == other.num_precursors and energy_deposition == other.energy_deposition and
         is_fissionable == other.is_fissionable and adjoint_mode == other.adjoint_mode;
}

void
XSProperties::MergeTransferMatricesProp(const XSProperties& other)
{
  if (transfer_matrix_props.size() != other.transfer_matrix_props.size())
    throw std::runtime_error("TransferMatrices size mismatch.");
  for (std::size_t ell = 0; ell < transfer_matrix_props.size(); ++ell)
  {
    transfer_matrix_props[ell] += other.transfer_matrix_props[ell];
  }
}

void
XSProperties::InsertPrecursor(const std::pair<double, std::vector<double>>& precursor)
{
  const auto is_identical = [&precursor](const auto& other)
  { return PrecursorsAreIdentical(precursor, other); };
  if (not std::ranges::any_of(precursors, is_identical))
    precursors.insert(precursor);
}

void
XSProperties::MergePrecursors(const XSProperties& other)
{
  for (const auto& precursor : other.precursors)
    InsertPrecursor(precursor);
}

Interpolator::Interpolator(const CartesianGrid& grid,
                           const std::vector<std::shared_ptr<MultiGroupXS>>& xs,
                           std::uint64_t flag)
  : grid_(grid), flag_(flag)
{
  if (grid.GetSize() != xs.size())
    throw std::runtime_error("Size mismatch.");

  // check XS properties
  if (xs.empty())
    return;
  if (not xs.front())
    throw std::runtime_error("Cross-section pointer cannot be null.");
  ref_xs_prop_ = XSProperties(*xs.front());
  for (std::size_t i = 1; i < xs.size(); ++i)
  {
    if (not xs[i])
      throw std::runtime_error("Cross-section pointer cannot be null.");
    XSProperties xs_prop(*xs[i]);
    if (xs_prop != ref_xs_prop_)
      throw std::runtime_error("Cross-section properties mismatch.");
    ref_xs_prop_.MergeTransferMatricesProp(xs_prop);
    ref_xs_prop_.MergePrecursors(xs_prop);
  }
  // compute itemsize and accumulate copying actions
  std::size_t itemsize = 0;
  std::vector<std::function<double*(double*, const MultiGroupXS&)>> getters;
  if (flag_ & XSType::Total)
  {
    itemsize += ref_xs_prop_.num_groups;
    getters.emplace_back(
      [](double* dest, const MultiGroupXS& xs)
      { return std::copy_n(xs.GetSigmaTotal().data(), xs.GetNumGroups(), dest); });
  }
  if (flag_ & XSType::Absorption)
  {
    itemsize += ref_xs_prop_.num_groups;
    getters.emplace_back(
      [](double* dest, const MultiGroupXS& xs)
      { return std::copy_n(xs.GetSigmaAbsorption().data(), xs.GetNumGroups(), dest); });
  }
  if (ref_xs_prop_.is_fissionable)
  {
    if (flag_ & XSType::Fission)
    {
      itemsize += ref_xs_prop_.num_groups;
      getters.emplace_back(
        [](double* dest, const MultiGroupXS& xs)
        { return std::copy_n(xs.GetSigmaFission().data(), xs.GetNumGroups(), dest); });
    }
    if (flag_ & XSType::NuFission)
    {
      itemsize += ref_xs_prop_.num_groups;
      getters.emplace_back(
        [](double* dest, const MultiGroupXS& xs)
        { return std::copy_n(xs.GetNuSigmaF().data(), xs.GetNumGroups(), dest); });
    }
    if (flag_ & XSType::Chi)
    {
      itemsize += ref_xs_prop_.num_groups;
      getters.emplace_back([](double* dest, const MultiGroupXS& xs)
                           { return std::copy_n(xs.GetChi().data(), xs.GetNumGroups(), dest); });
    }
    if (flag_ & XSType::ProductionMatrix)
    {
      itemsize += static_cast<std::size_t>(ref_xs_prop_.num_groups * ref_xs_prop_.num_groups);
      getters.emplace_back(
        [](double* dest, const MultiGroupXS& xs)
        {
          for (unsigned int g = 0; g < xs.GetNumGroups(); ++g)
            dest = std::copy_n(xs.GetProductionMatrix()[g].data(), xs.GetNumGroups(), dest);
          return dest;
        });
    }
    if (ref_xs_prop_.num_precursors > 0)
    {
      if (flag_ & XSType::NuPromptFission)
      {
        itemsize += ref_xs_prop_.num_groups;
        getters.emplace_back(
          [](double* dest, const MultiGroupXS& xs)
          { return std::copy_n(xs.GetNuPromptSigmaF().data(), xs.GetNumGroups(), dest); });
      }
      if (flag_ & XSType::NuDelayedFission)
      {
        itemsize += ref_xs_prop_.num_groups;
        getters.emplace_back(
          [](double* dest, const MultiGroupXS& xs)
          { return std::copy_n(xs.GetNuDelayedSigmaF().data(), xs.GetNumGroups(), dest); });
      }
      if (flag_ & XSType::Precursor)
      {
        itemsize += ref_xs_prop_.precursors.size();
        getters.emplace_back(
          [this](double* dest, const MultiGroupXS& xs)
          {
            for (const auto& precursor_prop : ref_xs_prop_.precursors)
              *(dest++) = GetFractionalYield(xs.GetPrecursors(), precursor_prop);
            return dest;
          });
      }
    }
  }

  if ((not ref_xs_prop_.transfer_matrix_props.empty()) and (flag_ & Transfer))
  {
    for (const auto& transfer_matrix_prop : ref_xs_prop_.transfer_matrix_props)
      itemsize += transfer_matrix_prop.present_indices.size();
    getters.emplace_back(
      [this](double* dest, const MultiGroupXS& xs)
      {
        for (std::size_t ell = 0; ell < ref_xs_prop_.transfer_matrix_props.size(); ++ell)
        {
          const auto& transfer_matrix_prop = ref_xs_prop_.transfer_matrix_props[ell];
          for (const auto& present_idx : transfer_matrix_prop.present_indices)
          {
            *dest = xs.GetTransferMatrices()[ell].GetValueIJ(present_idx.first, present_idx.second);
            dest++;
          }
        }
        return dest;
      });
  }

  // copy data from XS into contiguous buffer
  xs_data_ = RT_NDArray(grid.GetShape(), itemsize);
  std::array<std::uint32_t, MAX_DIM> nd_index_buffer{};
  std::span<std::uint32_t> nd_index(nd_index_buffer.data(), grid.GetShape().size());
  for (std::size_t i = 0; i < xs.size(); ++i)
  {
    const MultiGroupXS& xs_i = *(xs[i]);
    ContiguousToNDIndex(i, grid.GetShape(), nd_index);
    auto dest_view = xs_data_.Get(nd_index);

    double* dest = dest_view.data();
    for (const auto& getter : getters)
    {
      dest = getter(dest, xs_i);
    }
  }
}

MultiGroupXS
Interpolator::Evaluate(std::span<double> state_point)
{
  const auto shape = grid_.GetShape();
  if (state_point.size() != shape.size())
    throw std::invalid_argument("State point rank does not match the interpolation grid rank.");
  for (std::size_t d = 0; d < shape.size(); ++d)
  {
    const auto grid_d = grid_.GetGrid(static_cast<std::uint32_t>(d));
    const double state_value = state_point[d];
    if (state_value < grid_d.front() or state_value > grid_d.back())
      throw std::out_of_range("State point is outside the interpolation grid.");
  }

  MultiGroupXS result(ref_xs_prop_.num_groups,
                      ref_xs_prop_.scattering_order,
                      static_cast<unsigned int>(ref_xs_prop_.precursors.size()),
                      ref_xs_prop_.is_fissionable,
                      ref_xs_prop_.adjoint_mode);
  if (not ref_xs_prop_.energy_deposition.empty())
    result.GetEnergyDeposition() = ref_xs_prop_.energy_deposition;

  std::vector<double> contiguous_result = EvaluateContiguous(state_point);
  const double* data = contiguous_result.data();
  if (flag_ & XSType::Total)
  {
    result.GetSigmaTotal() = std::vector<double>(data, data + ref_xs_prop_.num_groups);
    data += ref_xs_prop_.num_groups;
  }
  if (flag_ & XSType::Absorption)
  {
    result.GetSigmaAbsorption() = std::vector<double>(data, data + ref_xs_prop_.num_groups);
    data += ref_xs_prop_.num_groups;
  }
  if (ref_xs_prop_.is_fissionable)
  {
    if (flag_ & XSType::Fission)
    {
      result.GetSigmaFission() = std::vector<double>(data, data + ref_xs_prop_.num_groups);
      data += ref_xs_prop_.num_groups;
    }
    if (flag_ & XSType::NuFission)
    {
      result.GetNuSigmaF() = std::vector<double>(data, data + ref_xs_prop_.num_groups);
      data += ref_xs_prop_.num_groups;
    }
    if (flag_ & XSType::Chi)
    {
      result.GetChi() = std::vector<double>(data, data + ref_xs_prop_.num_groups);
      data += ref_xs_prop_.num_groups;
    }
    if (flag_ & XSType::ProductionMatrix)
    {
      result.GetProductionMatrix().resize(ref_xs_prop_.num_groups);
      for (unsigned int g = 0; g < ref_xs_prop_.num_groups; ++g)
      {
        result.GetProductionMatrix()[g].assign(data, data + ref_xs_prop_.num_groups);
        data += ref_xs_prop_.num_groups;
      }
    }
    if (ref_xs_prop_.num_precursors > 0)
    {
      if (flag_ & XSType::NuPromptFission)
      {
        result.GetNuPromptSigmaF() = std::vector<double>(data, data + ref_xs_prop_.num_groups);
        data += ref_xs_prop_.num_groups;
      }
      if (flag_ & XSType::NuDelayedFission)
      {
        result.GetNuDelayedSigmaF() = std::vector<double>(data, data + ref_xs_prop_.num_groups);
        data += ref_xs_prop_.num_groups;
      }
      if (flag_ & XSType::Precursor)
      {
        result.GetPrecursors().reserve(ref_xs_prop_.precursors.size());
        for (const auto& precursor_prop : ref_xs_prop_.precursors)
          result.GetPrecursors().push_back(
            MultiGroupXS::Precursor(precursor_prop.first, *(data++), precursor_prop.second));
      }
    }
  }
  if ((not ref_xs_prop_.transfer_matrix_props.empty()) and (flag_ & Transfer))
  {
    result.GetTransferMatrices().reserve(ref_xs_prop_.transfer_matrix_props.size());
    for (std::size_t ell = 0; ell < ref_xs_prop_.transfer_matrix_props.size(); ++ell)
    {
      const auto& transfer_matrix_prop = ref_xs_prop_.transfer_matrix_props[ell];
      result.GetTransferMatrices().push_back(Reconstruct(transfer_matrix_prop, data));
    }
  }

  result.ComputeDiffusionParameters();

  return result;
}

} // namespace opensn
