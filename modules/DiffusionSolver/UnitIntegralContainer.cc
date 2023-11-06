#include "modules/DiffusionSolver/UnitIntegralContainer.h"

#include "framework/math/SpatialDiscretization/CellMappings/CellMapping.h"
#include "framework/math/SpatialDiscretization/FiniteElement/QuadraturePointData.h"

#include "framework/mesh/Cell/cell.h"

namespace chi_diffusion
{

UnitIntegralContainer::UnitIntegralContainer(MatDbl IntV_gradShapeI_gradShapeJ,
                                             MatVec3 IntV_shapeI_gradshapeJ,
                                             MatDbl IntV_shapeI_shapeJ,
                                             VecDbl IntV_shapeI,
                                             VecVec3 IntV_gradshapeI,
                                             std::vector<MatDbl> IntS_shapeI_shapeJ,
                                             std::vector<VecDbl> IntS_shapeI,
                                             std::vector<MatVec3> IntS_shapeI_gradshapeJ,
                                             std::vector<std::vector<int>> face_dof_mappings,
                                             size_t num_nodes)
  : IntV_gradShapeI_gradShapeJ_(std::move(IntV_gradShapeI_gradShapeJ)),
    IntV_shapeI_gradshapeJ_(std::move(IntV_shapeI_gradshapeJ)),
    IntV_shapeI_shapeJ_(std::move(IntV_shapeI_shapeJ)),
    IntV_shapeI_(std::move(IntV_shapeI)),
    IntV_gradshapeI_(std::move(IntV_gradshapeI)),
    IntS_shapeI_shapeJ_(std::move(IntS_shapeI_shapeJ)),
    IntS_shapeI_(std::move(IntS_shapeI)),
    IntS_shapeI_gradshapeJ_(std::move(IntS_shapeI_gradshapeJ)),
    face_dof_mappings_(std::move(face_dof_mappings)),
    num_nodes_(num_nodes)
{
}

UnitIntegralContainer
UnitIntegralContainer::Make(const chi_math::CellMapping& cell_mapping)
{
  typedef chi_math::finite_element::VolumetricQuadraturePointData VolQPData;
  typedef chi_math::finite_element::SurfaceQuadraturePointData FaceQPData;

  VolQPData internal_data = cell_mapping.MakeVolumetricQuadraturePointData();
  std::vector<FaceQPData> faces_qp_data;
  for (size_t f = 0; f < cell_mapping.ReferenceCell().faces_.size(); ++f)
    faces_qp_data.push_back(cell_mapping.MakeSurfaceQuadraturePointData(f));

  const auto n_dof_per_cell = internal_data.NumNodes();

  MatDbl IntV_gradshapeI_gradshapeJ(n_dof_per_cell, VecDbl(n_dof_per_cell));
  MatVec3 IntV_shapeI_gradshapeJ(n_dof_per_cell, VecVec3(n_dof_per_cell));
  MatDbl IntV_shapeI_shapeJ(n_dof_per_cell, VecDbl(n_dof_per_cell));
  VecDbl IntV_shapeI(n_dof_per_cell);
  VecVec3 IntV_gradshapeI(n_dof_per_cell);

  std::vector<MatDbl> IntS_shapeI_shapeJ(faces_qp_data.size());
  std::vector<VecDbl> IntS_shapeI(faces_qp_data.size());
  std::vector<MatVec3> IntS_shapeI_gradshapeJ(faces_qp_data.size());

  //  volume integrals
  for (unsigned int i = 0; i < n_dof_per_cell; ++i)
  {
    for (unsigned int j = 0; j < n_dof_per_cell; ++j)
    {
      for (const auto& qp : internal_data.QuadraturePointIndices())
      {
        IntV_gradshapeI_gradshapeJ[i][j] +=
          internal_data.ShapeGrad(i, qp).Dot(internal_data.ShapeGrad(j, qp)) *
          internal_data.JxW(qp);

        IntV_shapeI_gradshapeJ[i][j] +=
          internal_data.ShapeValue(i, qp) * internal_data.ShapeGrad(j, qp) * internal_data.JxW(qp);

        IntV_shapeI_shapeJ[i][j] +=
          internal_data.ShapeValue(i, qp) * internal_data.ShapeValue(j, qp) * internal_data.JxW(qp);
      } // for qp
    }   // for j

    for (const auto& qp : internal_data.QuadraturePointIndices())
    {
      IntV_shapeI[i] += internal_data.ShapeValue(i, qp) * internal_data.JxW(qp);

      IntV_gradshapeI[i] += internal_data.ShapeGrad(i, qp) * internal_data.JxW(qp);
    } // for qp
  }   // for i

  //  surface integrals
  for (size_t f = 0; f < faces_qp_data.size(); ++f)
  {
    IntS_shapeI_shapeJ[f].resize(n_dof_per_cell, VecDbl(n_dof_per_cell));
    IntS_shapeI[f].resize(n_dof_per_cell);
    IntS_shapeI_gradshapeJ[f].resize(n_dof_per_cell, VecVec3(n_dof_per_cell));

    for (unsigned int i = 0; i < n_dof_per_cell; ++i)
    {
      for (unsigned int j = 0; j < n_dof_per_cell; ++j)
      {
        for (const auto& qp : faces_qp_data[f].QuadraturePointIndices())
        {
          IntS_shapeI_shapeJ[f][i][j] += faces_qp_data[f].ShapeValue(i, qp) *
                                         faces_qp_data[f].ShapeValue(j, qp) *
                                         faces_qp_data[f].JxW(qp);

          IntS_shapeI_gradshapeJ[f][i][j] += faces_qp_data[f].ShapeValue(i, qp) *
                                             faces_qp_data[f].ShapeGrad(j, qp) *
                                             faces_qp_data[f].JxW(qp);
        } // for qp
      }   // for j

      for (const auto& qp : faces_qp_data[f].QuadraturePointIndices())
      {
        IntS_shapeI[f][i] += faces_qp_data[f].ShapeValue(i, qp) * faces_qp_data[f].JxW(qp);
      } // for qp
    }   // for i
  }     // for f

  //  unit integral data
  UnitIntegralContainer ui_data(IntV_gradshapeI_gradshapeJ,
                                IntV_shapeI_gradshapeJ,
                                IntV_shapeI_shapeJ,
                                IntV_shapeI,
                                IntV_gradshapeI,
                                IntS_shapeI_shapeJ,
                                IntS_shapeI,
                                IntS_shapeI_gradshapeJ,
                                cell_mapping.GetFaceNodeMappings(),
                                n_dof_per_cell);

  return ui_data;
}

double
UnitIntegralContainer::IntV_gradShapeI_gradShapeJ(unsigned int i, unsigned int j) const
{
  double value;
  auto& row_I = IntV_gradShapeI_gradShapeJ_.at(i);
  value = row_I.at(j);
  return value;
}
chi_mesh::Vector3
UnitIntegralContainer::IntV_shapeI_gradshapeJ(unsigned int i, unsigned int j) const
{
  chi_mesh::Vector3 value;
  auto& row_I = IntV_shapeI_gradshapeJ_.at(i);
  value = row_I.at(j);
  return value;
}
double
UnitIntegralContainer::IntV_shapeI_shapeJ(unsigned int i, unsigned int j) const
{
  double value;
  auto& row_I = IntV_shapeI_shapeJ_.at(i);
  value = row_I.at(j);
  return value;
}
double
UnitIntegralContainer::IntV_shapeI(unsigned int i) const
{
  double value = IntV_shapeI_.at(i);
  return value;
}
chi_mesh::Vector3
UnitIntegralContainer::IntV_gradshapeI(unsigned int i) const
{
  chi_mesh::Vector3 value;
  value = IntV_gradshapeI_.at(i);
  return value;
}
double
UnitIntegralContainer::IntS_shapeI_shapeJ(unsigned int face, unsigned int i, unsigned int j) const
{
  double value;
  auto& face_data = IntS_shapeI_shapeJ_.at(face);
  auto& rowI = face_data.at(i);
  value = rowI.at(j);
  return value;
}

double
UnitIntegralContainer::IntS_shapeI(unsigned int face, unsigned int i) const
{
  double value;
  auto& face_data = IntS_shapeI_.at(face);
  value = face_data.at(i);
  return value;
}

chi_mesh::Vector3
UnitIntegralContainer::IntS_shapeI_gradshapeJ(unsigned int face,
                                              unsigned int i,
                                              unsigned int j) const
{
  chi_mesh::Vector3 value;
  auto& face_data = IntS_shapeI_gradshapeJ_.at(face);
  auto& rowI = face_data.at(i);
  value = rowI.at(j);
  return value;
}

} // namespace chi_diffusion
