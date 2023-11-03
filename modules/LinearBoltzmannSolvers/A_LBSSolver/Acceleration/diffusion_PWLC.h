#pragma once

#include "diffusion.h"

namespace lbs::acceleration
{

class DiffusionPWLCSolver : public lbs::acceleration::DiffusionSolver
{
public:
  DiffusionPWLCSolver(std::string text_name,
                      const chi_math::SpatialDiscretization& sdm,
                      const chi_math::UnknownManager& uk_man,
                      std::map<uint64_t, BoundaryCondition> bcs,
                      MatID2XSMap map_mat_id_2_xs,
                      const std::vector<UnitCellMatrices>& unit_cell_matrices,
                      bool verbose);

  // 02c
  void AssembleAand_b(const std::vector<double>& q_vector) override;
  // 02d
  void Assemble_b(const std::vector<double>& q_vector) override;
  void Assemble_b(Vec petsc_q_vector) override;
};

} // namespace lbs::acceleration


