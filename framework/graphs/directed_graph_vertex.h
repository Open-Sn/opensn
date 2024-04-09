// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/graphs/graph.h"
#include <map>
#include <set>

namespace opensn
{

/**General implementation of a directed-graph vertex.*/
struct GraphVertex
{
  size_t id;
  void* context;

  std::set<size_t> us_edge;
  std::set<size_t> ds_edge;

  std::map<size_t, double> us_weights;
  std::map<size_t, double> ds_weights;

  GraphVertex(size_t id, void* context) : id(id), context(context) {}

  explicit GraphVertex(size_t id) : id(id), context(nullptr) {}

  GraphVertex(const GraphVertex& vertex)
  {
    this->id = vertex.id;
    this->context = vertex.context;

    us_edge = vertex.us_edge;
    ds_edge = vertex.ds_edge;
  }

  GraphVertex(GraphVertex&& vertex) noexcept
  {
    this->id = vertex.id;
    this->context = vertex.context;

    us_edge = vertex.us_edge;
    ds_edge = vertex.ds_edge;

    vertex.context = nullptr;
  }

  GraphVertex& operator=(const GraphVertex& vertex)
  {
    if (this == &vertex)
      return *this;

    this->id = vertex.id;
    this->context = vertex.context;

    us_edge = vertex.us_edge;
    ds_edge = vertex.ds_edge;

    return *this;
  }

  GraphVertex& operator=(GraphVertex&& vertex) noexcept
  {
    this->id = vertex.id;
    this->context = vertex.context;

    us_edge = vertex.us_edge;
    ds_edge = vertex.ds_edge;

    vertex.context = nullptr;

    return *this;
  }

  bool operator==(const GraphVertex& other) const { return other.id == this->id; }
};

} // namespace opensn
