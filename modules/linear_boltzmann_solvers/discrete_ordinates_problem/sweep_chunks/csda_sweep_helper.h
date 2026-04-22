// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include "framework/data_types/dense_matrix.h"
#include "framework/data_types/vector.h"
#include "framework/materials/multi_group_xs/multi_group_xs.h"
#include "framework/math/spatial_discretization/cell_mappings/cell_mapping.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/lbs_view.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/angle_set/angle_set.h"

#include <cmath>
#include <vector>

namespace opensn
{

struct CSDAMaterialData
{
  bool enabled = false;
  const std::vector<double>* stopping_power = nullptr;
  std::vector<double> delta_e;
};

inline constexpr double kCSDATolerance = 1.0e-12;

/// Context for one coupled local CSDA solve:
/// N nodal angular-flux unknowns plus one cellwise slope unknown.
struct CSDALocalSolveContext
{
  size_t gsg = 0;
  size_t gs_gi = 0;
  double sigma_tg = 0.0;
  double cell_volume = 0.0;
  const DenseMatrix<double>* Amat = nullptr;
  const DenseMatrix<double>* M = nullptr;
  const Vector<double>* v = nullptr;
  const Cell* cell = nullptr;
  const CellMapping* cell_mapping = nullptr;
  const std::vector<FaceOrientation>* face_orientations = nullptr;
  const std::vector<double>* face_mu_values = nullptr;
  const std::vector<Vector<double>>* int_f_shape_i = nullptr;
  const CellLBSView* cell_transport_view = nullptr;
  const AngleSet* angle_set = nullptr;
  int ni_preloc_face_counter = -1;
  const CSDAMaterialData* csda_data = nullptr;
  DenseMatrix<double>* Aext = nullptr;
  Vector<double>* b_ext = nullptr;
  std::vector<double>* psiE_gsg = nullptr;
};

inline CSDAMaterialData
MakeCSDAMaterialData(const MultiGroupXS& xs)
{
  CSDAMaterialData data;
  const auto& stopping_power = xs.GetStoppingPower();
  if (stopping_power.empty())
    return data;

  data.enabled = true;
  data.stopping_power = &stopping_power;
  data.delta_e = xs.GetDeltaE();
  return data;
}

inline bool
IsCSDAActive(const CSDAMaterialData& data, const size_t global_g)
{
  return data.enabled and data.stopping_power and global_g < data.stopping_power->size() and
         global_g < data.delta_e.size() and
         std::abs((*data.stopping_power)[global_g]) > kCSDATolerance;
}

template <class PsiEPtrAt>
inline void
StoreOutgoingCSDASlope(const bool should_store,
                       const std::vector<double>& psiE_gsg,
                       PsiEPtrAt psiE_ptr_at)
{
  if (not should_store)
    return;

  double* psiE = psiE_ptr_at();
  if (not psiE)
    return;

  for (size_t gsg = 0; gsg < psiE_gsg.size(); ++gsg)
    psiE[gsg] = psiE_gsg[gsg];
}

template <class SourceAt, class FluxAt, class FluxSet, class UpwindSlopeAt>
inline void
SolveCSDALocalSystem(const CSDALocalSolveContext& ctx,
                     SourceAt source_at,
                     FluxAt flux_at,
                     FluxSet flux_set,
                     UpwindSlopeAt upwind_slope_at)
{
  const auto& csda_data = *ctx.csda_data;
  auto& Aext = *ctx.Aext;
  auto& b_ext = *ctx.b_ext;
  auto& psiE_gsg = *ctx.psiE_gsg;
  const auto& Amat = *ctx.Amat;
  const auto& M = *ctx.M;
  const auto& v = *ctx.v;
  const auto& cell = *ctx.cell;
  const auto& cell_mapping = *ctx.cell_mapping;
  const auto& face_orientations = *ctx.face_orientations;
  const auto& face_mu_values = *ctx.face_mu_values;
  const auto& int_f_shape_i = *ctx.int_f_shape_i;
  const auto& cell_transport_view = *ctx.cell_transport_view;
  const size_t gsg = ctx.gsg;
  const size_t gs_gi = ctx.gs_gi;
  const double sigma_tg = ctx.sigma_tg;
  const double cell_volume = ctx.cell_volume;
  const int ni_preloc_face_counter = ctx.ni_preloc_face_counter;

  const double Sg = (*csda_data.stopping_power)[gs_gi + gsg];
  const double dEg = csda_data.delta_e[gs_gi + gsg];
  const size_t cell_num_nodes = cell_mapping.GetNumNodes();
  const size_t cell_num_faces = cell.faces.size();

  for (size_t i = 0; i < cell_num_nodes; ++i)
  {
    double temp = 0.0;
    for (size_t j = 0; j < cell_num_nodes; ++j)
    {
      const double Mij = M(i, j);
      Aext(i, j) = Amat(i, j) + Mij * (sigma_tg + Sg / dEg);
      temp += Mij * source_at(j);
    }

    double csda_gm1 = 0.0;
    if (gsg > 0)
    {
      const double Sgm1 = (*csda_data.stopping_power)[gs_gi + gsg - 1];
      const double dEgm1 = csda_data.delta_e[gs_gi + gsg - 1];
      double sumM = 0.0;
      for (size_t j = 0; j < cell_num_nodes; ++j)
        sumM += M(i, j) * flux_at(gsg - 1, j);
      csda_gm1 = (Sgm1 / dEgm1) * sumM - Sgm1 * v(i) * psiE_gsg[gsg - 1];
    }

    Aext(i, cell_num_nodes) = -Sg * v(i);
    b_ext(i) = flux_at(gsg, i) + temp + csda_gm1;
  }

  double slope_stream_coeff = 0.0;
  double slope_rhs_stream = 0.0;
  {
    int in_face_counter_slope = -1;
    int preloc_face_counter_slope = ni_preloc_face_counter;
    for (size_t f = 0; f < cell_num_faces; ++f)
    {
      const auto orientation = face_orientations[f];
      const bool is_incoming = (orientation == FaceOrientation::INCOMING);
      const bool is_outgoing = (orientation == FaceOrientation::OUTGOING);
      if (not is_incoming and not is_outgoing)
        continue;

      const auto& cell_face = cell.faces[f];
      const bool is_local_face = cell_transport_view.IsFaceLocal(f);
      const bool is_boundary_face = not cell_face.has_neighbor;

      if (is_incoming)
      {
        if (is_local_face)
          ++in_face_counter_slope;
        else if (not is_boundary_face)
          ++preloc_face_counter_slope;
      }

      const size_t num_face_nodes = cell_mapping.GetNumFaceNodes(f);
      double face_area = 0.0;
      for (size_t fi = 0; fi < num_face_nodes; ++fi)
      {
        const int i = cell_mapping.MapFaceNode(f, fi);
        face_area += int_f_shape_i[f](i);
      }

      if (is_outgoing)
      {
        slope_stream_coeff += face_mu_values[f] * face_area;
      }
      else
      {
        double psiE_upwind = 0.0;
        if (is_local_face)
        {
          psiE_upwind =
            upwind_slope_at(true, false, f, in_face_counter_slope, preloc_face_counter_slope, gsg);
        }
        else if (not is_boundary_face)
        {
          psiE_upwind =
            upwind_slope_at(false, false, f, in_face_counter_slope, preloc_face_counter_slope, gsg);
        }
        else if (ctx.angle_set->GetBoundaries().at(cell_face.neighbor_id)->IsReflecting())
        {
          psiE_upwind =
            upwind_slope_at(false, true, f, in_face_counter_slope, preloc_face_counter_slope, gsg);
        }
        slope_rhs_stream += (-face_mu_values[f]) * face_area * psiE_upwind;
      }
    }
  }

  const double coeff_row = 3.0 * Sg / (dEg * dEg);
  for (size_t j = 0; j < cell_num_nodes; ++j)
    Aext(cell_num_nodes, j) = v(j) * coeff_row;

  double rhs_csda_gm1 = 0.0;
  if (gsg > 0)
  {
    const double Sgm1 = (*csda_data.stopping_power)[gs_gi + gsg - 1];
    const double dEgm1 = csda_data.delta_e[gs_gi + gsg - 1];
    double sum_v_psi = 0.0;
    for (size_t j = 0; j < cell_num_nodes; ++j)
      sum_v_psi += v(j) * flux_at(gsg - 1, j);
    rhs_csda_gm1 =
      (3.0 / dEg) * Sgm1 * ((1.0 / dEgm1) * sum_v_psi - psiE_gsg[gsg - 1] * cell_volume);
  }

  Aext(cell_num_nodes, cell_num_nodes) =
    slope_stream_coeff + sigma_tg * cell_volume + (3.0 * Sg * cell_volume / dEg);
  b_ext(cell_num_nodes) = slope_rhs_stream + rhs_csda_gm1;

  GaussElimination(Aext, b_ext, cell_num_nodes + 1);

  for (size_t i = 0; i < cell_num_nodes; ++i)
    flux_set(gsg, i, b_ext(i));
  psiE_gsg[gsg] = b_ext(cell_num_nodes);
}

} // namespace opensn
