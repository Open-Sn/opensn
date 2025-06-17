// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "physics/problems/linear_boltzmann/response_evaluator/response_evaluator.h"
#include "physics/problems/linear_boltzmann/lbs_problem/point_source/point_source.h"
#include "physics/problems/linear_boltzmann/lbs_problem/volumetric_source/volumetric_source.h"
#include "physics/problems/linear_boltzmann/lbs_problem/io/lbs_problem_io.h"
#include "physics/solvers/solver.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "framework/logging/log.h"
#include "framework/object_factory.h"
#include "mpicpp-lite/mpicpp-lite.h"

namespace mpi = mpicpp_lite;

namespace opensn
{

OpenSnRegisterObjectInNamespace(lbs, ResponseEvaluator);

InputParameters
ResponseEvaluator::GetInputParameters()
{
  InputParameters params;
  params.SetGeneralDescription(
    "A utility class for evaluating responses using precomputed adjoint solutions "
    "and arbitrary forward sources.");
  params.SetDocGroup("LBSUtilities");

  params.AddRequiredParameter<std::shared_ptr<Problem>>("lbs_problem",
                                                        "A handle to an existing LBS problem.");
  params.AddOptionalParameterBlock(
    "options", ParameterBlock(), "The specification of adjoint buffers and forward to use.");
  params.LinkParameterToBlock("options", "response::OptionsBlock");

  return params;
}

std::shared_ptr<ResponseEvaluator>
ResponseEvaluator::Create(const ParameterBlock& params)
{
  auto& factory = opensn::ObjectFactory::GetInstance();
  return factory.Create<ResponseEvaluator>("lbs::ResponseEvaluator", params);
}

ResponseEvaluator::ResponseEvaluator(const InputParameters& params)
  : lbs_problem_(std::dynamic_pointer_cast<LBSProblem>(
      params.GetParamValue<std::shared_ptr<Problem>>("lbs_problem")))
{
  if (params.IsParameterValid("options"))
  {
    auto options = GetOptionsBlock();
    options.AssignParameters(params.GetParam("options"));
    SetOptions(options);
  }
}

InputParameters
ResponseEvaluator::GetOptionsBlock()
{
  InputParameters params;
  params.SetGeneralDescription("A block of options for the response evaluator for adding adjoint "
                               "buffers and defining forward sources.");
  params.SetDocGroup("LBSResponseEvaluator");

  params.AddOptionalParameterArray(
    "buffers", {}, "An array of tables containing adjoint buffer specifications.");
  params.LinkParameterToBlock("buffers", "response::BufferOptionsBlock");

  params.AddOptionalParameter("clear_sources", false, "A flag to clear existing sources.");
  params.AddOptionalParameterBlock(
    "sources", ParameterBlock(), "An array of tables containing source specification information.");
  params.LinkParameterToBlock("sources", "response::SourceOptionsBlock");

  return params;
}

void
ResponseEvaluator::SetOptions(const InputParameters& params)
{
  if (params.IsParameterValid("buffers"))
  {
    const auto& user_buffer_params = params.GetParam("buffers");
    user_buffer_params.RequireBlockTypeIs(ParameterBlockType::ARRAY);
    for (int p = 0; p < user_buffer_params.GetNumParameters(); ++p)
    {
      auto buffer_params = GetBufferOptionsBlock();
      buffer_params.AssignParameters(user_buffer_params.GetParam(p));
      SetBufferOptions(buffer_params);
    }
  }

  if (params.IsParameterValid("clear_sources"))
    if (params.GetParamValue<bool>("clear_sources"))
    {
      material_sources_.clear();
      point_sources_.clear();
      volumetric_sources_.clear();
      boundary_sources_.clear();
    }

  if (params.IsParameterValid("sources"))
  {
    auto source_params = GetSourceOptionsBlock();
    source_params.AssignParameters(params.GetParam("sources"));
    SetSourceOptions(source_params);
  }
}

InputParameters
ResponseEvaluator::GetBufferOptionsBlock()
{
  InputParameters params;
  params.SetGeneralDescription("Options for adding adjoint buffers to the response evaluator.");
  params.SetDocGroup("LBSResponseEvaluator");

  params.AddRequiredParameter<std::string>(
    "name",
    "A name given to the buffer to identify it when querying the response evaluation routine.");
  params.AddRequiredParameterBlock(
    "file_prefixes",
    "A table containing file prefixes for flux moments and angular flux binary files. "
    "These are keyed by \"flux_moments\" and \"angular_fluxes\", respectively.");

  return params;
}

void
ResponseEvaluator::SetBufferOptions(const InputParameters& input)
{
  auto params = opensn::ResponseEvaluator::GetBufferOptionsBlock();
  params.AssignParameters(input);

  const auto name = params.GetParamValue<std::string>("name");
  OpenSnInvalidArgumentIf(adjoint_buffers_.count(name) > 0,
                          "An adjoint buffer with name " + name + " already exists.");

  const auto prefixes = params.GetParam("file_prefixes");

  std::vector<double> phi;
  if (prefixes.Has("flux_moments"))
    LBSSolverIO::ReadFluxMoments(
      *lbs_problem_, prefixes.GetParamValue<std::string>("flux_moments"), false, phi);

  std::vector<std::vector<double>> psi;
  if (prefixes.Has("angular_fluxes"))
    LBSSolverIO::ReadAngularFluxes(
      *lbs_problem_, prefixes.GetParamValue<std::string>("angular_fluxes"), psi);

  adjoint_buffers_[name] = {phi, psi};
  log.Log0Verbose1() << "Adjoint buffer " << name << " added to the stack.";
}

InputParameters
ResponseEvaluator::GetSourceOptionsBlock()
{
  InputParameters params;
  params.SetGeneralDescription("A table of various forward source specifications.");
  params.SetDocGroup("LBSResponseEvaluator");

  params.AddOptionalParameterArray(
    "material", {}, "An array of tables containing material source specifications.");
  params.LinkParameterToBlock("material", "response::MaterialSourceOptionsBlock");

  params.AddOptionalParameterArray<std::shared_ptr<PointSource>>(
    "point", {}, "An array of tables containing point source handles.");

  params.AddOptionalParameterArray<std::shared_ptr<VolumetricSource>>(
    "volumetric", {}, "An array of tables containing volumetric source handles.");

  params.AddOptionalParameterArray(
    "boundary", {}, "An array of tables containing boundary source specifications.");
  params.LinkParameterToBlock("boundary", "BoundaryOptionsBlock");

  return params;
}

void
ResponseEvaluator::SetSourceOptions(const InputParameters& input)
{
  auto params = ResponseEvaluator::GetSourceOptionsBlock();
  params.AssignParameters(input);

  params.RequireBlockTypeIs(ParameterBlockType::BLOCK);

  // Add material sources
  if (params.Has("material"))
  {
    const auto& user_msrc_params = params.GetParam("material");
    for (int p = 0; p < user_msrc_params.GetNumParameters(); ++p)
    {
      auto msrc_params = GetMaterialSourceOptionsBlock();
      msrc_params.AssignParameters(user_msrc_params.GetParam(p));
      SetMaterialSourceOptions(msrc_params);
    }
  }

  // Add point sources
  if (params.Has("point"))
  {
    const auto& user_psrc_params = params.GetParam("point");
    for (int p = 0; p < user_psrc_params.GetNumParameters(); ++p)
    {
      point_sources_.push_back(
        user_psrc_params.GetParam(p).GetValue<std::shared_ptr<PointSource>>());
      point_sources_.back()->Initialize(*lbs_problem_);
    }
  }

  // Add volumetric sources
  if (params.Has("volumetric"))
  {
    const auto& user_dsrc_params = params.GetParam("volumetric");
    for (int p = 0; p < user_dsrc_params.GetNumParameters(); ++p)
    {
      volumetric_sources_.push_back(
        user_dsrc_params.GetParam(p).GetValue<std::shared_ptr<VolumetricSource>>());
      volumetric_sources_.back()->Initialize(*lbs_problem_);
    }
  }

  // Add boundary sources
  if (params.Has("boundary"))
  {
    const auto& user_bsrc_params = params.GetParam("boundary");
    for (int p = 0; p < user_bsrc_params.GetNumParameters(); ++p)
    {
      auto bsrc_params = LBSProblem::GetBoundaryOptionsBlock();
      bsrc_params.AssignParameters(user_bsrc_params.GetParam(p));
      SetBoundarySourceOptions(bsrc_params);
    }
  }
}

InputParameters
ResponseEvaluator::GetMaterialSourceOptionsBlock()
{
  InputParameters params;
  params.SetGeneralDescription(
    "Options for adding material-based forward sources to the response evaluator.");
  params.SetDocGroup("LBSResponseEvaluator");

  params.AddRequiredParameter<int>("block_id", "The block id the source belongs to.");
  params.AddRequiredParameterArray("strength", "The group-wise material source strength.");

  return params;
}

void
ResponseEvaluator::SetMaterialSourceOptions(const InputParameters& params)
{
  const auto blkid = params.GetParamValue<int>("block_id");
  OpenSnInvalidArgumentIf(material_sources_.count(blkid) > 0,
                          "A material source for block id " + std::to_string(blkid) +
                            " already exists.");

  const auto values = params.GetParamVectorValue<double>("strength");
  OpenSnInvalidArgumentIf(values.size() != lbs_problem_->GetNumGroups(),
                          "The number of material source values and groups "
                          "in the underlying solver do not match. "
                          "Expected " +
                            std::to_string(lbs_problem_->GetNumGroups()) + " but got " +
                            std::to_string(values.size()) + ".");

  material_sources_[blkid] = values;
  log.Log0Verbose1() << "Material source for block id " << blkid << " added to the stack.";
}

void
ResponseEvaluator::SetBoundarySourceOptions(const InputParameters& params)
{
  const auto bndry_name = params.GetParamValue<std::string>("name");
  const auto bndry_type = params.GetParamValue<std::string>("type");

  const auto bid = lbs_problem_->supported_boundary_names.at(bndry_name);
  if (bndry_type == "isotropic")
  {
    OpenSnInvalidArgumentIf(not params.Has("group_strength"),
                            "Parameter \"group_strength\" is required for "
                            "boundaries of type \"isotropic\".");
    params.RequireParameterBlockTypeIs("values", ParameterBlockType::ARRAY);

    boundary_sources_[bid] = {LBSBoundaryType::ISOTROPIC,
                              params.GetParamVectorValue<double>("group_strength")};
  }
  else
    log.Log0Warning() << "Unsupported boundary type. Skipping the entry.";
}

void
ResponseEvaluator::ClearForwardSources()
{
  material_sources_.clear();
  point_sources_.clear();
  volumetric_sources_.clear();
  boundary_sources_.clear();
}

void
ResponseEvaluator::AddResponseBuffers(const InputParameters& params)
{
  params.RequireBlockTypeIs(ParameterBlockType::ARRAY);
  for (size_t p = 0; p < params.GetNumParameters(); ++p)
  {
    auto spec = ResponseEvaluator::GetBufferOptionsBlock();
    spec.AssignParameters(params.GetParam(p));
    SetBufferOptions(spec);
  }
}

void
ResponseEvaluator::AddResponseSources(const InputParameters& params)
{
  auto spec = ResponseEvaluator::GetSourceOptionsBlock();
  spec.AssignParameters(params);
  SetSourceOptions(spec);
}

double
ResponseEvaluator::EvaluateResponse(const std::string& buffer) const
{
  const auto& buffer_data = adjoint_buffers_.at(buffer);
  const auto& phi_dagger = buffer_data.first;
  const auto& psi_dagger = buffer_data.second;

  OpenSnLogicalErrorIf(not material_sources_.empty() and phi_dagger.empty(),
                       "If material sources are present, adjoint flux moments "
                       "must be available for response evaluation.");
  OpenSnLogicalErrorIf(not point_sources_.empty() and phi_dagger.empty(),
                       "If point sources are set, adjoint flux moments "
                       "must be available for response evaluation.");
  OpenSnLogicalErrorIf(not volumetric_sources_.empty() and phi_dagger.empty(),
                       "if volumetric sources are set, adjoint flux moments "
                       "must be available for response evaluation.");
  OpenSnLogicalErrorIf(not boundary_sources_.empty() and psi_dagger.empty(),
                       "If boundary sources are set, adjoint angular fluxes "
                       "must be available for response evaluation.");

  const auto& grid = lbs_problem_->GetGrid();
  const auto& discretization = lbs_problem_->GetSpatialDiscretization();
  const auto& transport_views = lbs_problem_->GetCellTransportViews();
  const auto& unit_cell_matrices = lbs_problem_->GetUnitCellMatrices();
  const auto num_groups = lbs_problem_->GetNumGroups();

  double local_response = 0.0;

  // Material sources
  if (not material_sources_.empty())
  {
    for (const auto& cell : grid->local_cells)
    {
      const auto& cell_mapping = discretization.GetCellMapping(cell);
      const auto& transport_view = transport_views[cell.local_id];
      const auto& fe_values = unit_cell_matrices[cell.local_id];
      const auto num_cell_nodes = cell_mapping.GetNumNodes();

      if (material_sources_.count(cell.block_id) > 0)
      {
        const auto& src = material_sources_.at(cell.block_id);
        for (size_t i = 0; i < num_cell_nodes; ++i)
        {
          const auto dof_map = transport_view.MapDOF(i, 0, 0);
          const auto& V_i = fe_values.intV_shapeI(i);
          for (size_t g = 0; g < num_groups; ++g)
            local_response += src[g] * phi_dagger[dof_map + g] * V_i;
        }
      }
    } // for cell
  }   // if material sources

  // Boundary sources
  if (not boundary_sources_.empty())
  {
    size_t gs = 0;
    for (const auto& groupset : lbs_problem_->GetGroupsets())
    {
      const auto& uk_man = groupset.psi_uk_man_;
      const auto& quadrature = groupset.quadrature;
      const auto& num_gs_angles = quadrature->omegas.size();
      const auto& num_gs_groups = groupset.groups.size();

      for (const auto& cell : grid->local_cells)
      {
        const auto& cell_mapping = discretization.GetCellMapping(cell);
        const auto& fe_values = unit_cell_matrices[cell.local_id];

        size_t f = 0;
        for (const auto& face : cell.faces)
        {
          if (not face.has_neighbor and boundary_sources_.count(face.neighbor_id) > 0)
          {
            const auto bndry_id = face.neighbor_id;
            const auto num_face_nodes = cell_mapping.GetNumFaceNodes(f);
            for (size_t fi = 0; fi < num_face_nodes; ++fi)
            {
              const auto i = cell_mapping.MapFaceNode(f, fi);
              const auto& node = grid->vertices[cell.vertex_ids[i]];
              const auto& intF_shapeI = fe_values.intS_shapeI[f](i);

              const auto psi_bndry = EvaluateBoundaryCondition(bndry_id, node, groupset);

              for (size_t n = 0; n < num_gs_angles; ++n)
              {
                const auto& omega = quadrature->omegas[n];
                const auto mu = omega.Dot(face.normal);
                if (mu < 0.0)
                {
                  const auto& wt = quadrature->weights[n];
                  const auto weight = -mu * wt * intF_shapeI;
                  const auto dof_map = discretization.MapDOFLocal(cell, i, uk_man, n, 0);

                  for (size_t gsg = 0; gsg < num_gs_groups; ++gsg)
                    local_response +=
                      weight * psi_dagger[gs][dof_map + gsg] * psi_bndry[num_gs_groups * n + gsg];
                } // if outgoing
              }
            } // for face node fi
          }
          ++f;
        } // for face
      }   // for cell
      ++gs;
    } // for groupset
  }   // if boundary sources

  // Point sources
  for (const auto& point_source : point_sources_)
    for (const auto& subscriber : point_source->GetSubscribers())
    {
      const auto& cell = grid->local_cells[subscriber.cell_local_id];
      const auto& transport_view = transport_views[cell.local_id];

      const auto& src = point_source->GetStrength();
      const auto& vol_wt = subscriber.volume_weight;

      const auto num_cell_nodes = transport_view.GetNumNodes();
      for (size_t i = 0; i < num_cell_nodes; ++i)
      {
        const auto dof_map = transport_view.MapDOF(i, 0, 0);
        const auto& shape_val = subscriber.shape_values(i);
        for (size_t g = 0; g < num_groups; ++g)
          local_response += vol_wt * shape_val * src[g] * phi_dagger[dof_map + g];
      } // for node i
    }   // for subscriber

  // Volumetric sources
  for (const auto& volumetric_source : volumetric_sources_)
    for (const uint64_t local_id : volumetric_source->GetSubscribers())
    {
      const auto& cell = grid->local_cells[local_id];
      const auto& transport_view = transport_views[cell.local_id];
      const auto& fe_values = unit_cell_matrices[cell.local_id];
      const auto& nodes = discretization.GetCellNodeLocations(cell);

      const auto num_cell_nodes = transport_view.GetNumNodes();
      for (size_t i = 0; i < num_cell_nodes; ++i)
      {
        const auto& V_i = fe_values.intV_shapeI(i);
        const auto dof_map = transport_view.MapDOF(i, 0, 0);
        const auto& vals = (*volumetric_source)(cell, nodes[i], num_groups);
        for (size_t g = 0; g < num_groups; ++g)
          local_response += vals[g] * phi_dagger[dof_map + g] * V_i;
      }
    }

  double global_response = 0.0;
  mpi_comm.all_reduce(local_response, global_response, mpi::op::sum<double>());
  return global_response;
}

std::vector<double>
ResponseEvaluator::EvaluateBoundaryCondition(const uint64_t boundary_id,
                                             const Vector3& node,
                                             const LBSGroupset& groupset,
                                             const double) const
{
  const auto num_gs_angles = groupset.quadrature->omegas.size();
  const auto num_gs_groups = groupset.groups.size();
  const auto first_group = groupset.groups.front().id;

  std::vector<double> psi;
  const auto& bc = boundary_sources_.at(boundary_id);
  if (bc.type == LBSBoundaryType::ISOTROPIC)
  {
    for (size_t n = 0; n < num_gs_angles; ++n)
      for (size_t gsg = 0; gsg < num_gs_groups; ++gsg)
        psi.emplace_back(bc.isotropic_mg_source[first_group + gsg]);
    return psi;
  }
  OpenSnLogicalError("Unexpected behavior. Unsupported boundary condition encountered.");
}

} // namespace opensn
