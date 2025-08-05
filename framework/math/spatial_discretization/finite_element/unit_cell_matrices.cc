// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/spatial_discretization/finite_element/unit_cell_matrices.h"
#include "framework/math/spatial_discretization/spatial_discretization.h"
#include "framework/math/spatial_weight_function.h"
#include "framework/math/spatial_discretization/finite_element/finite_element_data.h"

namespace opensn
{

UnitCellMatrices
ComputeUnitCellIntegrals(const SpatialDiscretization& sdm,
                         const Cell& cell,
                         const CoordinateSystemType coord_sys)
{
  auto swf = SpatialWeightFunction::FromCoordinateType(coord_sys);
  const auto& cell_mapping = sdm.GetCellMapping(cell);
  const size_t cell_num_faces = cell.faces.size();
  const size_t cell_num_nodes = cell_mapping.GetNumNodes();
  const auto fe_vol_data = cell_mapping.MakeVolumetricFiniteElementData();

  DenseMatrix<double> IntV_gradshapeI_gradshapeJ(cell_num_nodes, cell_num_nodes, 0.0);
  DenseMatrix<Vector3> IntV_shapeI_gradshapeJ(cell_num_nodes, cell_num_nodes);
  DenseMatrix<double> IntV_shapeI_shapeJ(cell_num_nodes, cell_num_nodes, 0.0);
  Vector<double> IntV_shapeI(cell_num_nodes, 0.);
  std::vector<DenseMatrix<double>> IntS_shapeI_shapeJ(cell_num_faces);
  std::vector<DenseMatrix<Vector3>> IntS_shapeI_gradshapeJ(cell_num_faces);
  std::vector<Vector<double>> IntS_shapeI(cell_num_faces);

  // Volume integrals
  for (unsigned int i = 0; i < cell_num_nodes; ++i)
  {
    for (unsigned int j = 0; j < cell_num_nodes; ++j)
    {
      for (const auto& qp : fe_vol_data.GetQuadraturePointIndices())
      {
        IntV_gradshapeI_gradshapeJ(i, j) +=
          (*swf)(fe_vol_data.QPointXYZ(qp)) *
          fe_vol_data.ShapeGrad(i, qp).Dot(fe_vol_data.ShapeGrad(j, qp)) *
          fe_vol_data.JxW(qp); // K-matrix

        IntV_shapeI_gradshapeJ(i, j) +=
          (*swf)(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.ShapeValue(i, qp) *
          fe_vol_data.ShapeGrad(j, qp) * fe_vol_data.JxW(qp); // G-matrix

        IntV_shapeI_shapeJ(i, j) += (*swf)(fe_vol_data.QPointXYZ(qp)) *
                                    fe_vol_data.ShapeValue(i, qp) * fe_vol_data.ShapeValue(j, qp) *
                                    fe_vol_data.JxW(qp); // M-matrix
      } // for qp
    } // for j

    for (const auto& qp : fe_vol_data.GetQuadraturePointIndices())
    {
      IntV_shapeI(i) +=
        (*swf)(fe_vol_data.QPointXYZ(qp)) * fe_vol_data.ShapeValue(i, qp) * fe_vol_data.JxW(qp);
    } // for qp
  } // for i

  //  surface integrals
  for (size_t f = 0; f < cell_num_faces; ++f)
  {
    const auto fe_srf_data = cell_mapping.MakeSurfaceFiniteElementData(f);
    IntS_shapeI_shapeJ[f] = DenseMatrix<double>(cell_num_nodes, cell_num_nodes, 0.0);
    IntS_shapeI[f] = Vector<double>(cell_num_nodes, 0.);
    IntS_shapeI_gradshapeJ[f] = DenseMatrix<Vector3>(cell_num_nodes, cell_num_nodes);

    for (unsigned int i = 0; i < cell_num_nodes; ++i)
    {
      for (unsigned int j = 0; j < cell_num_nodes; ++j)
      {
        for (const auto& qp : fe_srf_data.GetQuadraturePointIndices())
        {
          IntS_shapeI_shapeJ[f](i, j) += (*swf)(fe_srf_data.QPointXYZ(qp)) *
                                         fe_srf_data.ShapeValue(i, qp) *
                                         fe_srf_data.ShapeValue(j, qp) * fe_srf_data.JxW(qp);
          IntS_shapeI_gradshapeJ[f](i, j) += (*swf)(fe_srf_data.QPointXYZ(qp)) *
                                             fe_srf_data.ShapeValue(i, qp) *
                                             fe_srf_data.ShapeGrad(j, qp) * fe_srf_data.JxW(qp);
        } // for qp
      } // for j

      for (const auto& qp : fe_srf_data.GetQuadraturePointIndices())
      {
        IntS_shapeI[f](i) +=
          (*swf)(fe_srf_data.QPointXYZ(qp)) * fe_srf_data.ShapeValue(i, qp) * fe_srf_data.JxW(qp);
      } // for qp
    } // for i
  } // for f

  return UnitCellMatrices{IntV_gradshapeI_gradshapeJ,
                          IntV_shapeI_gradshapeJ,
                          IntV_shapeI_shapeJ,
                          IntV_shapeI,

                          IntS_shapeI_shapeJ,
                          IntS_shapeI_gradshapeJ,
                          IntS_shapeI};
}

} // namespace opensn
