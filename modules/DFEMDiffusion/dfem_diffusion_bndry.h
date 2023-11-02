#ifndef DFEM_DIFFUSION_BOUNDARY_H
#define DFEM_DIFFUSION_BOUNDARY_H

#include <array>

namespace dfem_diffusion
{
  class Boundary;
  
  enum class BoundaryType : int
  {
    Reflecting = 1,
    Dirichlet  = 2,
    Neumann    = 3,
    Robin      = 4,
    Vacuum     = 5
  };
}

//###################################################################
/**Parent class for diffusion boundaries*/
class dfem_diffusion::Boundary
{
  public :
  BoundaryType type_ = BoundaryType::Dirichlet;

  std::array<double, 3> values_ = {0.,0.,0.};
};

#endif //DFEM_DIFFUSION_BOUNDARY_H