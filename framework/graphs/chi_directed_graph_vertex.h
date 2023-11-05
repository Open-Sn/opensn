#pragma once

#include "opensn/framework/graphs/chi_graph.h"

#include <map>

//###################################################################
/**General implementation of a directed-graph vertex.*/
struct chi::GraphVertex
{
  size_t id;
  void* context;

  std::set<size_t> us_edge;
  std::set<size_t> ds_edge;

  std::map<size_t, double> us_weights;
  std::map<size_t, double> ds_weights;

  GraphVertex(size_t in_id, void* in_context) : id(in_id), context(in_context) {}

  explicit GraphVertex(size_t in_id) : id(in_id), context(nullptr) {}

  GraphVertex(const GraphVertex& in_v)
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;
  }

  GraphVertex(GraphVertex&& in_v) noexcept
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;

    in_v.context = nullptr;
  }

  GraphVertex& operator=(const GraphVertex& in_v)
  {
    if (this == &in_v) return *this;

    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;

    return *this;
  }

  GraphVertex& operator=(GraphVertex&& in_v) noexcept
  {
    this->id = in_v.id;
    this->context = in_v.context;

    us_edge = in_v.us_edge;
    ds_edge = in_v.ds_edge;

    in_v.context = nullptr;

    return *this;
  }

  bool operator==(const GraphVertex& other) const { return other.id == this->id; }
};
