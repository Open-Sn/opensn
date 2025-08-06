// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/solvers/keigen_solver.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "framework/field_functions/field_function_grid_based.h"
#include "framework/runtime.h"

namespace opensn
{

InputParameters
KEigenSolver::GetInputParameters()
{
  InputParameters params = Solver::GetInputParameters();
  params.AddRequiredParameter<std::shared_ptr<Problem>>("problem", "An existing lbs problem");
  return params;
}

KEigenSolver::KEigenSolver(const InputParameters& params)
  : Solver(params),
    do_problem_(params.GetSharedPtrParam<Problem, DiscreteOrdinatesProblem>("problem"))
{
}

void
KEigenSolver::InitializePowerFieldFunction()
{
  auto& options = do_problem_->GetOptions();
  std::shared_ptr<SpatialDiscretization> sdm(&do_problem_->GetSpatialDiscretization());
  // Initialize power generation field function
  std::string prefix;
  if (options.field_function_prefix_option == "prefix")
  {
    prefix = options.field_function_prefix;
    if (not prefix.empty())
      prefix += "_";
  }
  if (options.field_function_prefix_option == "solver_name")
    prefix = GetName() + "_";

  power_field_func_ = std::make_shared<FieldFunctionGridBased>(
    prefix + "power_generation", sdm, Unknown(UnknownType::SCALAR));
}

std::shared_ptr<FieldFunctionGridBased>
KEigenSolver::GetPowerFieldFunction() const
{
  return power_field_func_;
}

void
KEigenSolver::UpdatePowerFieldFunctions()
{
  // Update power generation
  const auto local_node_count_ = do_problem_->GetLocalNodeCount();
  const auto& grid_ = do_problem_->GetGrid();
  const auto& sdm = do_problem_->GetSpatialDiscretization();
  const auto& unit_cell_matrices_ = do_problem_->GetUnitCellMatrices();
  const auto& block_id_to_xs_map_ = do_problem_->GetMatID2XSMap();
  const auto& groups_ = do_problem_->GetGroups();
  auto& phi_new_local_ = do_problem_->GetPhiNewLocal();
  const auto& options_ = do_problem_->GetOptions();
  const auto& phi_uk_man = do_problem_->GetUnknownManager();

  std::vector<double> data_vector_local(local_node_count_, 0.0);

  double local_total_power = 0.0;
  for (const auto& cell : grid_->local_cells)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t num_nodes = cell_mapping.GetNumNodes();

    const auto& Vi = unit_cell_matrices_[cell.local_id].intV_shapeI;

    const auto& xs = block_id_to_xs_map_.at(cell.block_id);

    if (not xs->IsFissionable())
      continue;

    for (size_t i = 0; i < num_nodes; ++i)
    {
      const int64_t imapA = sdm.MapDOFLocal(cell, i);
      const int64_t imapB = sdm.MapDOFLocal(cell, i, phi_uk_man, 0, 0);

      double nodal_power = 0.0;
      for (size_t g = 0; g < groups_.size(); ++g)
      {
        const double sigma_fg = xs->GetSigmaFission()[g];
        // const double kappa_g = xs->Kappa()[g];
        const double kappa_g = options_.power_default_kappa;

        nodal_power += kappa_g * sigma_fg * phi_new_local_[imapB + g];
      } // for g

      data_vector_local[imapA] = nodal_power;
      local_total_power += nodal_power * Vi(i);
    } // for node
  } // for cell

  if (options_.power_normalization > 0.0)
  {
    double global_total_power;
    mpi_comm.all_reduce(local_total_power, global_total_power, mpi::op::sum<double>());

    Scale(data_vector_local, options_.power_normalization / global_total_power);
  }

  // const size_t ff_index = power_gen_fieldfunc_local_handle_;

  // auto& ff_ptr = field_functions_.at(ff_index);
  power_field_func_->UpdateFieldVector(data_vector_local);
}

} // namespace opensn
