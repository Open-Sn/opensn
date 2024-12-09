// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/mesh/cell/cell.h"

namespace opensn
{

/// Stores references to global cells to enable an iterator
class LocalCellHandler
{
public:
  std::vector<std::unique_ptr<Cell>>& local_cells;

private:
  /// Constructor.
  explicit LocalCellHandler(std::vector<std::unique_ptr<Cell>>& cells) : local_cells(cells) {}

public:
  static LocalCellHandler Create(std::vector<std::unique_ptr<Cell>>& native_cells)
  {
    return LocalCellHandler(native_cells);
  }

  /// Returns a reference to a local cell, given a local cell index.
  Cell& operator[](uint64_t cell_local_index);

  /// Returns a const reference to a local cell, given a local cell index.
  const Cell& operator[](uint64_t cell_local_index) const;

  /**
   * Get the total number of local cells.
   * @return Size of the native cell vector.
   */
  size_t LocalCellCount() const { return local_cells.size(); }

  /// Internal iterator class.
  class Iterator
  {
  private:
    LocalCellHandler& handler_;
    size_t index_;

  public:
    Iterator(LocalCellHandler& handler, size_t index) : handler_(handler), index_(index) {}

    Iterator& operator++()
    {
      ++index_;
      return *this;
    }

    Cell& operator*() { return *(handler_.local_cells[index_]); }
    bool operator==(const Iterator& other) const { return index_ == other.index_; }
    bool operator!=(const Iterator& other) const { return index_ != other.index_; }
  };

  /// Internal const iterator class.
  class ConstIterator
  {
  private:
    const LocalCellHandler& handler_;
    size_t index_;

  public:
    ConstIterator(const LocalCellHandler& handler, size_t index) : handler_(handler), index_(index)
    {
    }

    ConstIterator operator++()
    {
      ++index_;
      return *this;
    }

    const Cell& operator*() { return *(handler_.local_cells[index_]); }
    bool operator==(const ConstIterator& other) const { return index_ == other.index_; }
    bool operator!=(const ConstIterator& other) const { return index_ != other.index_; }
  };

  Iterator begin() { return {*this, 0}; }

  Iterator end() { return {*this, local_cells.size()}; }

  ConstIterator begin() const { return {*this, 0}; }

  ConstIterator end() const { return {*this, local_cells.size()}; }

  size_t size() const { return local_cells.size(); }
};

} // namespace opensn
