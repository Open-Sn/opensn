// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/diffusion/boundary_closure.h"
#include "framework/mesh/mesh/cell.h"
#include "framework/mesh/mesh/mesh.h"
#include "framework/math/spatial_discretization/cell_mappings/cell_mapping.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"
#include "framework/math/spatial_weight_function.h"
#include "framework/utils/error.h"
#include <petscblaslapack.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace opensn
{

double
QuadratureHalfRangeMoment(const std::vector<Vector3>& omegas,
                          const std::vector<double>& weights,
                          const Vector3& normal)
{
  OpenSnInvalidArgumentIf(omegas.size() != weights.size() or omegas.empty(),
                          "QuadratureHalfRangeMoment: directions and weights must be non-empty "
                          "and of equal size.");
  double sum_w = 0.0;
  double sum_w_mu = 0.0;
  for (size_t d = 0; d < omegas.size(); ++d)
  {
    sum_w += weights[d];
    sum_w_mu += weights[d] * std::fabs(omegas[d].Dot(normal));
  }
  OpenSnInvalidArgumentIf(sum_w <= 0.0, "QuadratureHalfRangeMoment: weights must sum to > 0.");
  return 0.5 * sum_w_mu / sum_w;
}

double
FaceOrthogonalLength(const Cell& cell, unsigned int f)
{
  const double volume = cell.volume;
  const double face_area = cell.faces.at(f).area;
  const size_t num_faces = cell.faces.size();
  double surface_area = 0.0;
  for (const auto& face : cell.faces)
    surface_area += face.area;

  switch (cell.GetType())
  {
    case CellType::SLAB:
      return volume;
    case CellType::POLYGON:
      if (num_faces == 3)
        return 2.0 * volume / face_area;
      if (num_faces == 4)
        return volume / face_area;
      return 4.0 * volume / surface_area;
    case CellType::POLYHEDRON:
      if (num_faces == 4)
        return 3.0 * volume / face_area;
      if (num_faces == 6 and cell.vertex_ids.size() == 8)
        return volume / face_area;
      return 6.0 * volume / surface_area;
    default:
      OpenSnLogicalError("FaceOrthogonalLength: unsupported cell type.");
  }
}

double
FaceSideLength(const Mesh& grid, const Cell& cell, unsigned int f)
{
  const auto& face = cell.faces.at(f);
  const Vector3& xc = cell.centroid;
  switch (cell.GetType())
  {
    case CellType::SLAB:
      return cell.volume;
    case CellType::POLYGON:
    {
      // Side triangle (v0, v1, xc): height over the edge = 2 |T| / |e|.
      const Vector3& v0 = grid.GlobalVertex(face.vertex_ids.at(0));
      const Vector3& v1 = grid.GlobalVertex(face.vertex_ids.at(1));
      const double edge = (v1 - v0).Norm();
      const double twice_area = std::fabs((v0 - xc).Cross(v1 - xc).z);
      OpenSnLogicalErrorIf(edge <= 0.0 or twice_area <= 0.0,
                           "FaceSideLength: degenerate face or side triangle (cell centroid on the "
                           "line of a face).");
      return 2.0 * twice_area / edge;
    }
    case CellType::POLYHEDRON:
    {
      // Side tets (v_k, v_k+1, face centroid, cell centroid): height over their face triangle
      // = 3 |S_k| / |T_k|.
      const Vector3& xf = face.centroid;
      const size_t nv = face.vertex_ids.size();
      double hmin = std::numeric_limits<double>::max();
      for (size_t k = 0; k < nv; ++k)
      {
        const Vector3& a = grid.GlobalVertex(face.vertex_ids[k]);
        const Vector3& b = grid.GlobalVertex(face.vertex_ids[(k + 1) % nv]);
        const Vector3 normal_x2 = (a - xf).Cross(b - xf);
        const double twice_tri = normal_x2.Norm();
        if (twice_tri <= 0.0)
          continue;
        // 6 |S_k| = |normal_x2 . (xc - xf)|, 2 |T_k| = twice_tri
        hmin = std::min(hmin, std::fabs(normal_x2.Dot(xc - xf)) / twice_tri);
      }
      OpenSnLogicalErrorIf(hmin == std::numeric_limits<double>::max() or hmin <= 0.0,
                           "FaceSideLength: degenerate face or side tetrahedron (cell centroid on "
                           "the plane of a face).");
      return 2.0 * hmin;
    }
    default:
      OpenSnLogicalError("FaceSideLength: unsupported cell type.");
  }
}

double
MaxEigenvalueOnStiffnessRange(const DenseMatrix<double>& N, const DenseMatrix<double>& K)
{
  const size_t n = K.Rows();
  OpenSnInvalidArgumentIf(
    K.Columns() != n or N.Rows() != n or N.Columns() != n,
    "MaxEigenvalueOnStiffnessRange: matrices must be square and of equal size.");

  // K is singular on the constants and N vanishes on them, so the maximum of v'Nv / v'Kv over the
  // range of K equals the maximum generalized eigenvalue of N against K + alpha 1 1' (any alpha >
  // 0).
  double trace = 0.0;
  for (size_t i = 0; i < n; ++i)
    trace += K(i, i);
  const double alpha = trace / static_cast<double>(n * n);

  std::vector<double> a(n * n), b(n * n);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
    {
      a[j * n + i] = 0.5 * (N(i, j) + N(j, i));
      b[j * n + i] = 0.5 * (K(i, j) + K(j, i)) + alpha;
    }

  const PetscBLASInt itype = 1;
  const char jobz = 'N';
  const char uplo = 'U';
  const auto n_blas = static_cast<PetscBLASInt>(n);
  const PetscBLASInt lwork = std::max(PetscBLASInt{1}, 3 * n_blas);
  std::vector<double> w(n), work(lwork);
  PetscBLASInt info = 0;
  LAPACKsygv_(&itype,
              &jobz,
              &uplo,
              &n_blas,
              a.data(),
              &n_blas,
              b.data(),
              &n_blas,
              w.data(),
              work.data(),
              &lwork,
              &info);
  OpenSnLogicalErrorIf(
    info != 0,
    "MaxEigenvalueOnStiffnessRange: LAPACK sygv failed (info = " + std::to_string(info) +
      "); the stiffness matrix may have a null space beyond the constants.");
  return *std::max_element(w.begin(), w.end());
}

double
FaceTraceConstant(const CellMapping& cell_mapping,
                  const Cell& cell,
                  unsigned int f,
                  const DenseMatrix<double>& stiffness,
                  CoordinateSystemType coord_sys)
{
  const auto swf = SpatialWeightFunction::FromCoordinateType(coord_sys);
  const auto num_nodes = cell_mapping.GetNumNodes();
  const auto& normal = cell.faces.at(f).normal;
  const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);
  DenseMatrix<double> N(num_nodes, num_nodes, 0.0);
  for (size_t i = 0; i < num_nodes; ++i)
    for (size_t j = 0; j < num_nodes; ++j)
      for (size_t qp : fe_srf_data.GetQuadraturePointIndices())
        N(i, j) += (*swf)(fe_srf_data.QPointXYZ(qp)) * normal.Dot(fe_srf_data.ShapeGrad(i, qp)) *
                   normal.Dot(fe_srf_data.ShapeGrad(j, qp)) * fe_srf_data.JxW(qp);
  return FaceOrthogonalLength(cell, f) * MaxEigenvalueOnStiffnessRange(N, stiffness);
}

double
ClosureRTheta(
  double beta, double sigma_t_h, double c_grad, unsigned int num_closure_faces, double s0)
{
  // Smooth switch from the Marshak closure (theta = 1) on thin cells to the half-range form
  // (theta = 0) on thick ones.
  double theta = 1.0;
  if (sigma_t_h > 0.0)
    theta = std::min(1.0, std::pow(s0 / sigma_t_h, 2));

  // Coercivity floor: the root in [0, 1] of (1 - t)^2 / (1 + t) = r.
  if (c_grad > 0.0 and num_closure_faces > 0)
  {
    const double r = 6.0 * beta * sigma_t_h / (c_grad * static_cast<double>(num_closure_faces));
    if (r < 1.0)
    {
      const double theta_coer = 0.5 * ((2.0 + r) - std::sqrt(r * r + 8.0 * r));
      theta = std::max(theta, std::clamp(theta_coer, 0.0, 1.0));
    }
  }
  return std::clamp(theta, 0.0, 1.0);
}

} // namespace opensn
