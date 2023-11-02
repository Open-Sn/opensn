#include "diffusion.h"

#include "math/SpatialDiscretization/SpatialDiscretization.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

namespace lbs::acceleration
{

// ###################################################################
/**Default constructor.*/
DiffusionSolver::DiffusionSolver(
  std::string text_name,
  const chi_math::SpatialDiscretization& sdm,
  const chi_math::UnknownManager& uk_man,
  std::map<uint64_t, BoundaryCondition> bcs,
  MatID2XSMap map_mat_id_2_xs,
  const std::vector<UnitCellMatrices>& unit_cell_matrices,
  const bool verbose,
  const bool requires_ghosts)
  : text_name_(std::move(text_name)),
    grid_(sdm.Grid()),
    sdm_(sdm),
    uk_man_(uk_man),
    bcs_(std::move(bcs)),
    mat_id_2_xs_map_(std::move(map_mat_id_2_xs)),
    unit_cell_matrices_(unit_cell_matrices),
    num_local_dofs_(static_cast<int64_t>(sdm_.GetNumLocalDOFs(uk_man_))),
    num_global_dofs_(static_cast<int64_t>(sdm_.GetNumGlobalDOFs(uk_man_))),
    A_(nullptr),
    rhs_(nullptr),
    ksp_(nullptr),
    requires_ghosts_(requires_ghosts)
{
  options.verbose = verbose;
}

// ###################################################################
/**Default destructor.*/
DiffusionSolver::~DiffusionSolver()
{
  MatDestroy(&A_);
  VecDestroy(&rhs_);
  KSPDestroy(&ksp_);
}

// ###################################################################
/**Returns the assigned text name.*/
std::string DiffusionSolver::TextName() const { return text_name_; }

// ###################################################################
/**Returns the right-hand side petsc vector.*/
const Vec& DiffusionSolver::RHS() const { return rhs_; }

// ###################################################################
/**Returns the assigned unknown structure.*/
const chi_math::UnknownManager& DiffusionSolver::UnknownStructure() const
{
  return uk_man_;
}

// ###################################################################
/**Returns the associated spatial discretization.*/
const chi_math::SpatialDiscretization&
DiffusionSolver::SpatialDiscretization() const
{
  return sdm_;
}

std::pair<size_t, size_t>
lbs::acceleration::DiffusionSolver::GetNumPhiIterativeUnknowns()
{
  return {sdm_.GetNumLocalDOFs(uk_man_), sdm_.GetNumGlobalDOFs(uk_man_)};
}

// ##################################################################
/**Adds to the right-hand side without applying spatial discretization.*/
void lbs::acceleration::DiffusionSolver::AddToRHS(
  const std::vector<double>& values)
{
  typedef unsigned int uint;
  typedef const int64_t cint64_t;
  const size_t num_local_dofs = sdm_.GetNumLocalDOFs(uk_man_);

  ChiInvalidArgumentIf(num_local_dofs != values.size(),
                       "Vector size mismatched with spatial discretization");

  const size_t num_unknowns = uk_man_.NumberOfUnknowns();

  for (const auto& cell : grid_.local_cells)
  {
    const auto& cell_mapping = sdm_.GetCellMapping(cell);

    for (size_t i = 0; i < cell_mapping.NumNodes(); ++i)
    {
      for (size_t u = 0; u < num_unknowns; ++u)
      {
        for (uint c = 0; c < uk_man_.GetUnknown(u).NumComponents(); ++c)
        {
          cint64_t dof_map_local = sdm_.MapDOFLocal(cell, i, uk_man_, u, c);
          cint64_t dof_map = sdm_.MapDOF(cell, i, uk_man_, u, c);

          VecSetValue(rhs_, dof_map, values[dof_map_local], ADD_VALUES);
        } // for component c
      }   // for unknown u
    }     // for node i
  }       // for cell

  VecAssemblyBegin(rhs_);
  VecAssemblyEnd(rhs_);
}

} // namespace lbs::acceleration
