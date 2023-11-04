#pragma once

#include "../chi_mesh.h"

namespace chi_mesh
{

enum class SurfaceMesherType
{
  Predefined = 1
};

/**
 * Base class for surface meshers.
 */
class SurfaceMesher
{
protected:
  const SurfaceMesherType type_;
  std::vector<double> xcuts_;
  std::vector<double> ycuts_;

public:
  SurfaceMesherType GetType() const { return type_; }
  void AddXCut(double x_cut) { xcuts_.push_back(x_cut); }
  void AddYCut(double y_cut) { ycuts_.push_back(y_cut); }
  explicit SurfaceMesher(SurfaceMesherType in_type) : type_(in_type) {}

  virtual void Execute();

  virtual ~SurfaceMesher() = default;
};

} // namespace chi_mesh
