#include "lbs_solver.h"

#include "math/SpatialDiscretization/FiniteElement/PiecewiseLinear/PieceWiseLinearDiscontinuous.h"
#include "mesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_runtime.h"
#include "chi_log.h"

#include "console/chi_console.h"

#include <iomanip>

void lbs::LBSSolver::InitializeSpatialDiscretization()
{
  using namespace chi_math::finite_element;
  Chi::log.Log() << "Initializing spatial discretization.\n";
  discretization_ = chi_math::spatial_discretization::PieceWiseLinearDiscontinuous::New(*grid_ptr_);

  ComputeUnitIntegrals();
}

void lbs::LBSSolver::ComputeUnitIntegrals()
{
  Chi::log.Log() << "Computing unit integrals.\n";
  const auto& sdm = *discretization_;

  //======================================== Define spatial weighting functions
  struct SpatialWeightFunction //SWF
  {
    virtual double operator()(const chi_mesh::Vector3& pt) const
    { return 1.0; }
    virtual ~SpatialWeightFunction() = default;
  };

  struct SphericalSWF : public SpatialWeightFunction
  {
    double operator()(const chi_mesh::Vector3& pt) const override
    { return pt[2]*pt[2]; }
  };

  struct CylindricalSWF : public SpatialWeightFunction
  {
    double operator()(const chi_mesh::Vector3& pt) const override
    { return pt[0]; }
  };

  auto swf_ptr = std::make_shared<SpatialWeightFunction>();
  if (options_.geometry_type == lbs::GeometryType::ONED_SPHERICAL)
    swf_ptr = std::make_shared<SphericalSWF>();
  if (options_.geometry_type == lbs::GeometryType::TWOD_CYLINDRICAL)
    swf_ptr = std::make_shared<CylindricalSWF>();

  auto ComputeCellUnitIntegrals = [&sdm](const chi_mesh::Cell& cell,
                                         const SpatialWeightFunction& swf)
  {
    const auto& cell_mapping = sdm.GetCellMapping(cell);
    const size_t cell_num_faces = cell.faces_.size();
    const size_t cell_num_nodes = cell_mapping.NumNodes();
    const auto vol_qp_data = cell_mapping.MakeVolumetricQuadraturePointData();

    MatDbl  IntV_gradshapeI_gradshapeJ(cell_num_nodes, VecDbl(cell_num_nodes));
    MatVec3 IntV_shapeI_gradshapeJ(cell_num_nodes, VecVec3(cell_num_nodes));
    MatDbl  IntV_shapeI_shapeJ(cell_num_nodes, VecDbl(cell_num_nodes));
    VecDbl  IntV_shapeI(cell_num_nodes);

    std::vector<MatDbl>  IntS_shapeI_shapeJ(cell_num_faces);
    std::vector<MatVec3> IntS_shapeI_gradshapeJ(cell_num_faces);
    std::vector<VecDbl>  IntS_shapeI(cell_num_faces);

    //Volume integrals
    for (unsigned int i = 0; i < cell_num_nodes; ++i)
    {
      for (unsigned int j = 0; j < cell_num_nodes; ++j)
      {
        for (const auto& qp : vol_qp_data.QuadraturePointIndices())
        {
          IntV_gradshapeI_gradshapeJ[i][j]
            += swf(vol_qp_data.QPointXYZ(qp)) *
               vol_qp_data.ShapeGrad(i, qp).Dot(vol_qp_data.ShapeGrad(j, qp)) *
               vol_qp_data.JxW(qp);  //K-matrix

          IntV_shapeI_gradshapeJ[i][j]
            += swf(vol_qp_data.QPointXYZ(qp)) *
               vol_qp_data.ShapeValue(i, qp) *
               vol_qp_data.ShapeGrad(j, qp) *
               vol_qp_data.JxW(qp);  //G-matrix

          IntV_shapeI_shapeJ[i][j]
            += swf(vol_qp_data.QPointXYZ(qp)) *
               vol_qp_data.ShapeValue(i, qp) *
               vol_qp_data.ShapeValue(j, qp) *
               vol_qp_data.JxW(qp);  //M-matrix
        }// for qp
      }// for j

      for (const auto& qp : vol_qp_data.QuadraturePointIndices())
      {
        IntV_shapeI[i]
          += swf(vol_qp_data.QPointXYZ(qp)) *
             vol_qp_data.ShapeValue(i, qp) * vol_qp_data.JxW(qp);
      }// for qp
    }//for i

    //  surface integrals
    for (size_t f = 0; f < cell_num_faces; ++f)
    {
      const auto faces_qp_data = cell_mapping.MakeSurfaceQuadraturePointData(f);
      IntS_shapeI_shapeJ[f].resize(cell_num_nodes, VecDbl(cell_num_nodes));
      IntS_shapeI[f].resize(cell_num_nodes);
      IntS_shapeI_gradshapeJ[f].resize(cell_num_nodes, VecVec3(cell_num_nodes));

      for (unsigned int i = 0; i < cell_num_nodes; ++i)
      {
        for (unsigned int j = 0; j < cell_num_nodes; ++j)
        {
          for (const auto& qp : faces_qp_data.QuadraturePointIndices())
          {
            IntS_shapeI_shapeJ[f][i][j]
              += swf(faces_qp_data.QPointXYZ(qp)) *
                 faces_qp_data.ShapeValue(i, qp) *
                 faces_qp_data.ShapeValue(j, qp) *
                 faces_qp_data.JxW(qp);
            IntS_shapeI_gradshapeJ[f][i][j]
              += swf(faces_qp_data.QPointXYZ(qp)) *
                 faces_qp_data.ShapeValue(i, qp) *
                 faces_qp_data.ShapeGrad(j, qp) *
                 faces_qp_data.JxW(qp);
          }// for qp
        }//for j

        for (const auto& qp : faces_qp_data.QuadraturePointIndices())
        {
          IntS_shapeI[f][i]
            += swf(faces_qp_data.QPointXYZ(qp)) *
               faces_qp_data.ShapeValue(i, qp) * faces_qp_data.JxW(qp);
        }// for qp
      }//for i
    }//for f

    return
      UnitCellMatrices{IntV_gradshapeI_gradshapeJ, //K-matrix
                       IntV_shapeI_gradshapeJ,     //G-matrix
                       IntV_shapeI_shapeJ,         //M-matrix
                       IntV_shapeI,                //Vi-vectors

                       IntS_shapeI_shapeJ,         //face M-matrices
                       IntS_shapeI_gradshapeJ,     //face G-matrices
                       IntS_shapeI};               //face Si-vectors
  };

  const size_t num_local_cells = grid_ptr_->local_cells.size();
  unit_cell_matrices_.resize(num_local_cells);

  for (const auto& cell : grid_ptr_->local_cells)
    unit_cell_matrices_[cell.local_id_] = ComputeCellUnitIntegrals(cell, *swf_ptr);

  const auto ghost_ids = grid_ptr_->cells.GetGhostGlobalIDs();
  for (uint64_t ghost_id : ghost_ids)
    unit_ghost_cell_matrices_[ghost_id] =
      ComputeCellUnitIntegrals(grid_ptr_->cells[ghost_id],*swf_ptr);

  //============================================= Assessing global unit cell
  //                                              matrix storage
  std::array<size_t,2> num_local_ucms = {unit_cell_matrices_.size(),
                                         unit_ghost_cell_matrices_.size()};
  std::array<size_t,2> num_globl_ucms = {0,0};

  MPI_Allreduce(num_local_ucms.data(), //sendbuf
                num_globl_ucms.data(), //recvbuf
                2, MPIU_SIZE_T,        //count+datatype
                MPI_SUM,               //operation
                Chi::mpi.comm);       //comm



  Chi::mpi.Barrier();
  Chi::log.Log()
  << "Ghost cell unit cell-matrix ratio: "
  << (double)num_globl_ucms[1]*100/(double)num_globl_ucms[0]
  << "%";
  Chi::log.Log()
    << "Cell matrices computed.                   Process memory = "
    << std::setprecision(3)
    << chi::Console::GetMemoryUsageInMB() << " MB";
}