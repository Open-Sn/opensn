#pragma once

#include <vector>
#include "opensn/framework/mesh/chi_mesh.h"

namespace chi_math
{
class CellMapping;
}

namespace chi_diffusion
{

class UnitIntegralContainer
{
public:
  typedef std::vector<double> VecDbl;
  typedef std::vector<VecDbl> MatDbl;
  typedef std::vector<chi_mesh::Vector3> VecVec3;
  typedef std::vector<VecVec3> MatVec3;

private:
  MatDbl IntV_gradShapeI_gradShapeJ_;
  MatVec3 IntV_shapeI_gradshapeJ_;
  MatDbl IntV_shapeI_shapeJ_;
  VecDbl IntV_shapeI_;
  VecVec3 IntV_gradshapeI_;

  std::vector<MatDbl> IntS_shapeI_shapeJ_;
  std::vector<VecDbl> IntS_shapeI_;
  std::vector<MatVec3> IntS_shapeI_gradshapeJ_;

  std::vector<std::vector<int>> face_dof_mappings_;
  size_t num_nodes_ = 0;

public:
  UnitIntegralContainer(MatDbl IntV_gradShapeI_gradShapeJ,
                        MatVec3 IntV_shapeI_gradshapeJ,
                        MatDbl IntV_shapeI_shapeJ,
                        VecDbl IntV_shapeI,
                        VecVec3 IntV_gradshapeI,
                        std::vector<MatDbl> IntS_shapeI_shapeJ,
                        std::vector<VecDbl> IntS_shapeI,
                        std::vector<MatVec3> IntS_shapeI_gradshapeJ,
                        std::vector<std::vector<int>> face_dof_mappings,
                        size_t num_nodes);

  static UnitIntegralContainer Make(const chi_math::CellMapping& cell_mapping);

  double IntV_gradShapeI_gradShapeJ(unsigned int i, unsigned int j) const;

  chi_mesh::Vector3 IntV_shapeI_gradshapeJ(unsigned int i, unsigned int j) const;
  double IntV_shapeI_shapeJ(unsigned int i, unsigned int j) const;

  double IntV_shapeI(unsigned int i) const;

  chi_mesh::Vector3 IntV_gradshapeI(unsigned int i) const;

  double IntS_shapeI_shapeJ(unsigned int face, unsigned int i, unsigned int j) const;

  double IntS_shapeI(unsigned int face, unsigned int i) const;

  chi_mesh::Vector3 IntS_shapeI_gradshapeJ(unsigned int face, unsigned int i, unsigned int j) const;

  int FaceDofMapping(size_t face, size_t face_node_index) const
  {
    auto& face_data = face_dof_mappings_.at(face);
    return face_data.at(face_node_index);
  }

  size_t NumNodes() const { return num_nodes_; }

  const MatDbl& GetIntV_gradShapeI_gradShapeJ() const { return IntV_gradShapeI_gradShapeJ_; }
  const MatVec3& GetIntV_shapeI_gradshapeJ() const { return IntV_shapeI_gradshapeJ_; }
  const MatDbl& GetIntV_shapeI_shapeJ() const { return IntV_shapeI_shapeJ_; }
  const VecDbl& GetIntV_shapeI() const { return IntV_shapeI_; }
  const VecVec3& GetIntV_gradshapeI() const { return IntV_gradshapeI_; }

  const std::vector<MatDbl>& GetIntS_shapeI_shapeJ() const { return IntS_shapeI_shapeJ_; }
  const std::vector<VecDbl>& GetIntS_shapeI() const { return IntS_shapeI_; }
  const std::vector<MatVec3>& GetIntS_shapeI_gradshapeJ() const { return IntS_shapeI_gradshapeJ_; }
};

} // namespace chi_diffusion
