#pragma once

#include "framework/graphs/directed_graph_vertex.h"
#include <stack>

namespace opensn
{

/**Simple implementation of a directed graph. This implementation was
 * considered to serve more versatile strategies with regards to grid
 * parallel partitioning.*/
class DirectedGraph
{
public:
  /**Allows semi-sane access to vertices even if
   * they are removed from the graph.*/
  class VertexAccessor
  {
  private:
    std::vector<GraphVertex> vertices_;
    std::vector<bool> vertex_valid_flags_;

  public:
    /** Adds a vertex to the graph with a supplied id.*/
    void AddVertex(size_t id, void* context);
    /** Adds a vertex to the graph where the ID is assigned to
     * the number of vertices already loaded on the graph.
     * For example, if there are 3 vertices on the graph (with
     * IDs 0 through 2) then the next vertex (this one) will
     * be assigned and ID of 3.*/
    void AddVertex(void* context);
    /** Removes a vertex from the graph.*/
    void RemoveVertex(size_t v);

    /** Accesses a vertex from the graph.*/
    GraphVertex& operator[](size_t v);

    /**Internal iterator class for vertex accessor.*/
    class iterator
    {
    public:
      VertexAccessor& ref_block;
      size_t ref_element;

      iterator(VertexAccessor& in_block, size_t i) : ref_block(in_block), ref_element(i) {}

      iterator operator++()
      {
        iterator i = *this;
        bool stop = false;
        do
        {
          ++ref_element;
          if (ref_element >= ref_block.vertices_.size())
            stop = true;
          if (not stop)
            if (ref_block.vertex_valid_flags_[ref_element])
              stop = true;
        } while (not stop);
        return i;
      }
      iterator operator++(int junk)
      {
        bool stop = false;
        do
        {
          ++ref_element;
          if (ref_element >= ref_block.vertices_.size())
            stop = true;
          if (not stop)
            if (ref_block.vertex_valid_flags_[ref_element])
              stop = true;
        } while (not stop);
        return *this;
      }
      GraphVertex& operator*() { return ref_block.vertices_[ref_element]; }
      GraphVertex* operator->() { return &(ref_block.vertices_[ref_element]); }
      bool operator==(const iterator& rhs) const { return ref_element == rhs.ref_element; }
      bool operator!=(const iterator& rhs) const { return ref_element != rhs.ref_element; }
    };

    iterator begin()
    {
      size_t count = 0;
      if (vertex_valid_flags_.empty())
        return {*this, count};
      if (vertex_valid_flags_[count])
        return {*this, count};

      // Get next valid or end
      while (true)
      {
        if (count >= vertices_.size())
          return {*this, count};
        if (vertex_valid_flags_[count])
          return {*this, count};
        ++count;
      }
    }

    iterator end() { return {*this, vertices_.size()}; }

    size_t size() { return vertices_.size(); }

    size_t GetNumValid()
    {
      size_t count = 0;
      for (bool val : vertex_valid_flags_)
        if (val)
          ++count;

      return count;
    }

    void clear()
    {
      vertices_.clear();
      vertex_valid_flags_.clear();
    }
  };

  VertexAccessor vertices;

  /** Adds a vertex to the graph. By default <I>context</I> is
   * assumed to be nullptr.*/
  void AddVertex(size_t id, void* context = nullptr);
  /** Adds a vertex to the graph. By default <I>context</I> is
   * assumed to be nullptr and <I>id</I> is assumed to be assigned
   * automatically. In
   * the latter case the vertex id will be the same as the order
   * in which it was added (0,1,2,3,etc ... will have id's 0,1,2,3,etc)*/
  void AddVertex(void* context = nullptr);
  /** Removes a vertex from the graph. This method does not
   * free any context related data.*/
  void RemoveVertex(size_t v);
  /** Adds an edge to the graph. Range checks are supplied by the
   * vertex accessor.*/
  bool AddEdge(size_t from, size_t to, double weight = 1.0);
  /**Remove an edge from the graph. Range checks are supplied by the
   * vertex accessor.*/
  void RemoveEdge(size_t from, size_t to);

  size_t GetNumSinks()
  {
    size_t count = 0;
    for (auto& v : vertices)
      if (v.ds_edge.empty() and not v.us_edge.empty())
        ++count;
    return count;
  }

  size_t GetNumSources()
  {
    size_t count = 0;
    for (auto& v : vertices)
      if (v.us_edge.empty() and not v.ds_edge.empty())
        ++count;
    return count;
  }

private:
  /** Depth-First-Search main recursive algorithm. This is the recursive
   * portion of the method below this one
   * (opensn::DirectedGraph::DepthFirstSearch).*/
  void DFSAlgorithm(std::vector<size_t>& traversal, std::vector<bool>& visited, size_t cur_vid);

  /**SCC main recursive algorithm. This is the recursive call for the
   * method defined below this one
   * (opensn::DirectedGraph::FindStronglyConnectedConnections).*/
  void SCCAlgorithm(size_t u,
                    int& time,
                    std::vector<int>& disc,
                    std::vector<int>& low,
                    std::vector<bool>& on_stack,
                    std::stack<size_t>& stack,
                    std::vector<std::vector<size_t>>& SCCs);

public:
  /**Find strongly connected components. This method is the implementation
   * of Tarjan's algorithm [1].
   *
   * [1] Tarjan R.E. "Depth-first search and linear graph algorithms",
   *     SIAM Journal on Computing, 1972.
   *
   * It returns collections of vertices that form strongly connected
   * components excluding singletons.*/
  std::vector<std::vector<size_t>> FindStronglyConnectedComponents();

  /** Generates a topological sort. This method is the implementation
   * of Kahn's algorithm [1].
   *
   * [1] Kahn, Arthur B. (1962), "Topological sorting of large networks",
   *     Communications of the ACM, 5 (11): 558–562
   *
   * \return Returns the vertex ids sorted topologically. If this
   *         vector is empty the algorithm failed because it detected
   *         cyclic dependencies.*/
  std::vector<size_t> GenerateTopologicalSort();

  /**Finds a sequence that minimizes the Feedback Arc Set (FAS). This
   * algorithm implements the algorithm depicted in [1].
   *
   * [1] Eades P., Lin X., Smyth W.F., "Fast & Effective heuristic for
   *     the feedback arc set problem", Information Processing Letters,
   *     Volume 47. 1993.*/
  std::vector<size_t> FindApproxMinimumFAS();

  /**Prints the graph in Graphviz format.*/
  void PrintGraphviz(int location_mask = 0);

  /**Prints a sub-graph in Graphviz format.*/
  void PrintSubGraphviz(const std::vector<int>& verts_to_print, int location_mask = 0);

  std::vector<std::pair<size_t, size_t>> RemoveCyclicDependencies();

  /**Clears all the data structures associated with the graph.*/
  void Clear();

  ~DirectedGraph();
};
} // namespace opensn
