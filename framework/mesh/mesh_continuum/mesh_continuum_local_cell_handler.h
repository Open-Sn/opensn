// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/cell/cell.h"

namespace opensn
{

/// Stores references to global cells to enable an iterator
class LocalCellHandler
{
  friend class MeshContinuum;

public:
  std::vector<std::unique_ptr<Cell>>& native_cells;

private:
  /// Constructor.
  explicit LocalCellHandler(std::vector<std::unique_ptr<Cell>>& native_cells)
    : native_cells(native_cells)
  {
  }

public:
  /// Returns a reference to a local cell, given a local cell index.
  Cell& operator[](uint64_t cell_local_index);

  /// Returns a const reference to a local cell, given a local cell index.
  const Cell& operator[](uint64_t cell_local_index) const;

  /// Internal iterator class.
  class Iterator
  {
  public:
    LocalCellHandler& ref_block;
    size_t ref_element;

    Iterator(LocalCellHandler& block, size_t i) : ref_block(block), ref_element(i) {}

    Iterator operator++()
    {
      Iterator i = *this;
      ref_element++;
      return i;
    }
    Iterator operator++(int)
    {
      ref_element++;
      return *this;
    }

    Cell& operator*() { return *(ref_block.native_cells[ref_element]); }
    bool operator==(const Iterator& rhs) const { return ref_element == rhs.ref_element; }
    bool operator!=(const Iterator& rhs) const { return ref_element != rhs.ref_element; }
  };

  /// Internal const iterator class.
  class ConstIterator
  {
  public:
    const LocalCellHandler& ref_block;
    size_t ref_element;

    ConstIterator(const LocalCellHandler& block, size_t i) : ref_block(block), ref_element(i) {}

    ConstIterator operator++()
    {
      ConstIterator i = *this;
      ref_element++;
      return i;
    }
    ConstIterator operator++(int)
    {
      ref_element++;
      return *this;
    }

    const Cell& operator*() { return *(ref_block.native_cells[ref_element]); }
    bool operator==(const ConstIterator& rhs) const { return ref_element == rhs.ref_element; }
    bool operator!=(const ConstIterator& rhs) const { return ref_element != rhs.ref_element; }
  };

  Iterator begin() { return {*this, 0}; }

  Iterator end() { return {*this, native_cells.size()}; }

  ConstIterator begin() const { return {*this, 0}; }

  ConstIterator end() const { return {*this, native_cells.size()}; }

  size_t size() const { return native_cells.size(); }
};

} // namespace opensn
