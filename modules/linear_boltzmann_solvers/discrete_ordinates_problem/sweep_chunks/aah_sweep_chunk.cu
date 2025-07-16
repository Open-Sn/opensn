// SPDX-FileCopyrightText: 2025 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep_chunks/aah_sweep_chunk.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/discrete_ordinates_problem.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aahd_fluds.h"
#include "modules/linear_boltzmann_solvers/discrete_ordinates_problem/sweep/fluds/aah_fluds.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/memory_pinner.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/mesh_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/outflow_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/carrier/quadrature_carrier.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/storage.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/mesh_view.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/quadrature_view.h"
#include "modules/linear_boltzmann_solvers/lbs_problem/device/view/xs_view.h"
#include "framework/mesh/mesh_continuum/mesh_continuum.h"
#include "caliper/cali.h"
#include "caribou/caribou.h"

namespace crb = caribou;

namespace opensn
{

inline constexpr std::uint32_t matrix_size = max_dof * max_dof;

/// @brief Arguments for the sweep kernel.
struct SweepKernelArgs
{
  // mesh and quadrature
  const char* mesh_data;
  const char* quad_data;
  // source moments and phi
  const double* src_moment;
  double* phi;
  // angle set
  const std::uint32_t* directions;
  std::uint32_t angleset_size;
  // group set
  std::uint32_t num_groups;
  std::uint32_t groupset_start;
  std::uint32_t groupset_size;
  // fluds
  char* flud_data;
  std::uint32_t flud_stride;
  // level
  const std::uint32_t* level;
  std::uint32_t level_size;
  // batch
  std::uint32_t batch_size;
};

/// @brief Index for each thread.
struct Index
{
  /// @brief Constructor
  __device__ inline Index() {}

  __device__ inline void Compute(std::uint32_t thread_idx,
                                 const std::uint32_t& angleset_size,
                                 const std::uint32_t& groupset_size)
  {
    cell_idx = thread_idx / (groupset_size * angleset_size);
    group_idx = (thread_idx / angleset_size) % groupset_size;
    angle_idx = thread_idx % angleset_size;
  }

  /// @brief Index of the cell associated to the current thread in the current level vector.
  std::uint32_t cell_idx;
  /// @brief Index of the group associated to the current thread in the current groupset.
  std::uint32_t group_idx;
  /// @brief Index of the angle associated to the current thread in the current angleset.
  std::uint32_t angle_idx;
};

/// @brief Get the pointer to the face data of a given cell in the flatten FLUDs
__device__ inline std::int64_t*
GetFaceOffsetData(char* flud_data, int cell_level_idx)
{
  std::uint64_t* cell_offset = reinterpret_cast<std::uint64_t*>(flud_data);
  return reinterpret_cast<std::int64_t*>(flud_data + cell_offset[cell_level_idx]);
}

/**
 * @details Compute the source contribution for the angular flux and the contribution to the sweep
 * matrix from the transfer and the mass matrix.
 */
__device__ inline void
ComputeGMS(std::array<double, matrix_size>& sweep_matrix,
           std::array<double, max_dof>& psi,
           const std::uint32_t& cell_num_nodes,
           DirectionView& direction,
           CellView& cell,
           const double* src_moment,
           const std::uint32_t& groupset_start,
           const std::uint32_t& group_idx,
           const std::uint32_t& num_groups,
           const std::uint32_t& num_moments)
{
  // get sigmaT
  double sigma_t = cell.total_xs[groupset_start + group_idx];
  // compute source term
  std::array<double, max_dof> s;
  s.fill(0.0);
  src_moment += cell.phi_address + groupset_start + group_idx;
  for (std::uint32_t i = 0; i < cell_num_nodes; ++i)
  {
    double src_per_moment = 0.0;
    for (std::uint32_t m = 0; m < num_moments; ++m)
    {
      src_per_moment += direction.m2d[m] * (*src_moment);
      src_moment += num_groups;
    }
    s[i] = src_per_moment;
  }
  // add source, transfer and mass contribution
  double* A = sweep_matrix.data();
  const std::array<double, 4>* GM_data =
    reinterpret_cast<const std::array<double, 4>*>(cell.GM_data);
  for (std::uint32_t i = 0; i < cell_num_nodes; ++i)
  {
    for (std::uint32_t j = 0; j < cell_num_nodes; ++j)
    {
      std::array<double, 4> GM = *(GM_data++);
      // compute A += G * Omega + M * sigma_t
      A[j] += direction.omega[0] * GM[0] + direction.omega[1] * GM[1] + direction.omega[2] * GM[2] +
              sigma_t * GM[3];
      // compute psi += M @ s
      psi[i] += GM[3] * s[j];
    }
    A += max_dof;
  }
}

/// @brief Compute contribution from surface integral of the upwind angular flux.
__device__ inline void
ComputeSurfaceIntegral(std::array<double, matrix_size>& sweep_matrix,
                       std::array<double, max_dof>& psi,
                       const std::uint32_t& cell_num_nodes,
                       DirectionView& direction,
                       CellView& cell,
                       const char* flud_data,
                       const std::int64_t* face_offset_data,
                       const std::uint32_t& flud_stride)
{
  // loop over each face
  for (std::uint32_t f = 0; f < cell.num_faces; ++f)
  {
    // skip if not incoming face
    std::int64_t face_offset = face_offset_data[f];
    if (face_offset >= 0)
    {
      continue;
    }
    FaceView face;
    cell.GetFaceView(face, f);
    // compute scalar product with normal vector
    double mu = direction.omega[0] * face.normal[0] + direction.omega[1] * face.normal[1] +
                direction.omega[2] * face.normal[2];
    // get psi of the current angle and group, psi_s of different face node are separated by
    // flud_stride
    const double* upwind_psi = reinterpret_cast<const double*>(flud_data - face_offset);
    // compute surface integral
    for (std::uint32_t fi = 0; fi < face.num_face_nodes; ++fi)
    {
      std::uint32_t i = face.cell_mapping_data[fi];
      double* Ai = sweep_matrix.data() + i * max_dof;
      for (std::uint32_t fj = 0; fj < face.num_face_nodes; ++fj)
      {
        std::uint32_t j = face.cell_mapping_data[fj];
        double mu_Nij = -mu * face.M_surf_data[fi * face.num_face_nodes + fj];
        Ai[j] += mu_Nij;
        psi[i] += upwind_psi[fj * flud_stride] * mu_Nij;
      }
    }
  }
}

/// @brief Solve the sweep linear system using Gaussian elimination.
__device__ inline void
GaussianElimination(std::array<double, matrix_size>& sweep_matrix,
                    std::array<double, max_dof>& psi,
                    const std::uint32_t& cell_num_nodes)
{
  // forward elimination
  double* A_i = sweep_matrix.data();
  for (std::uint32_t i = 0; i < cell_num_nodes; ++i)
  {
    double inv_diag = 1.0 / A_i[i];
    // normalize the pivot row
    for (std::uint32_t j = i; j < cell_num_nodes; ++j)
    {
      A_i[j] *= inv_diag;
    }
    psi[i] *= inv_diag;
    // eliminate rows below
    double* A_k = A_i + max_dof;
    for (std::uint32_t k = i + 1; k < cell_num_nodes; ++k)
    {
      double factor = -A_k[i];
      for (std::uint32_t j = i; j < cell_num_nodes; ++j)
      {
        A_k[j] += factor * A_i[j];
      }
      psi[k] += factor * psi[i];
      A_k += max_dof;
    }
    A_i += max_dof;
  }
  // back substitution — row-wise access
  for (std::int32_t j = cell_num_nodes - 2; j >= 0; --j)
  {
    double* A_j = sweep_matrix.data() + j * max_dof;
    for (std::int32_t i = j + 1; i < cell_num_nodes; ++i)
    {
      psi[j] -= A_j[i] * psi[i];
    }
  }
}

/// @brief Record angular flux to downwind and compute outflow for boundary faces.
__device__ inline void
PutPsiToFludsAndOutflow(const std::array<double, max_dof>& psi,
                        DirectionView& direction,
                        CellView& cell,
                        char* flud_data,
                        const std::int64_t* face_offset_data,
                        const std::uint32_t& flud_stride,
                        const std::uint32_t& group,
                        const std::uint32_t& cell_local_id)
{
  // loop over each face
  for (std::uint32_t f = 0; f < cell.num_faces; ++f)
  {
    // skip if not outgoing face
    std::int64_t face_offset = face_offset_data[f];
    if (face_offset <= 0)
    {
      continue;
    }
    FaceView face;
    cell.GetFaceView(face, f);
    // get psi of the current angle and group, psi_s of different face node are separated by
    // flud_stride
    double* downwind_psi = reinterpret_cast<double*>(flud_data + face_offset);
    // loop over each face node
    for (std::uint32_t fi = 0; fi < face.num_face_nodes; ++fi)
    {
      std::uint32_t i = face.cell_mapping_data[fi];
      // put copy psi to FLUD
      downwind_psi[fi * flud_stride] = psi[i];
      // compute ouflow for boundary face
      if (face.outflow != nullptr)
      {
        double mu = direction.omega[0] * face.normal[0] + direction.omega[1] * face.normal[1] +
                    direction.omega[2] * face.normal[2];
        double outflow = direction.weight * mu * face.IntS_shapeI_data[fi] * psi[i];
        crb::atomic_add(face.outflow + group, outflow);
      }
    }
  }
}

/// @brief Compute the scalar flux.
__device__ void
ComputePhi(const std::array<double, max_dof>& psi,
           const std::uint32_t& cell_num_nodes,
           DirectionView& direction,
           CellView& cell,
           double* phi,
           const std::uint32_t& groupset_start,
           const std::uint32_t& group_idx,
           const std::uint32_t& num_groups,
           const std::uint32_t& num_moments)
{
  phi += cell.phi_address + groupset_start + group_idx;
  for (std::uint32_t i = 0; i < cell_num_nodes; ++i)
  {
    for (std::uint32_t m = 0; m < num_moments; ++m)
    {
      crb::atomic_add(phi, direction.d2m[m] * psi[i]);
      phi += num_groups;
    }
  }
}

/// @brief Save angular flux
__device__ inline void
RecordPsi(const std::array<double, max_dof>& psi,
          const std::uint32_t& cell_num_nodes,
          double* psi_device,
          Index& idx,
          const std::uint32_t& angleset_size,
          const std::uint32_t& groupset_size)
{
  // compute the stride
  std::uint32_t stride =
    idx.cell_idx * angleset_size * groupset_size + idx.angle_idx * groupset_size + idx.group_idx;
  stride *= max_dof;
  psi_device += stride;
  // set 0 for empty fluxes
  for (std::uint32_t i = 0; i < max_dof; ++i)
  {
    psi_device[i] = 0.0;
  }
  // record the flux
  for (std::uint32_t i = 0; i < cell_num_nodes; ++i)
  {
    psi_device[i] = psi[i];
  }
}

/// @brief Main sweep kernel.
__global__ void
SweepKernel(SweepKernelArgs args, double* psi_device)
{
  // compute index (cell, group, angle) from thread flatten index
  Index idx;
  {
    std::uint32_t thread_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_idx >= args.batch_size)
    {
      return;
    }
    idx.Compute(thread_idx, args.angleset_size, args.groupset_size);
  }
  // get corresponding cell
  std::uint32_t cell_local_idx = args.level[idx.cell_idx];
  CellView cell;
  MeshView(args.mesh_data).GetCellView(cell, cell_local_idx);
  // get corresponding direction and number of moments
  std::uint32_t num_moments;
  std::uint32_t direction_num = args.directions[idx.angle_idx];
  DirectionView direction;
  {
    QuadratureView quadrature(args.quad_data);
    num_moments = quadrature.num_moments;
    quadrature.GetDirectionView(direction, direction_num);
  }
  // compute face offset of the current cell in the FLUD
  const std::int64_t* face_offset_data = GetFaceOffsetData(args.flud_data, idx.cell_idx);
  char* flud_data =
    args.flud_data + (idx.angle_idx * args.groupset_size + idx.group_idx) * sizeof(double);
  // initialize psi
  std::array<double, max_dof> psi;
  psi.fill(0.0);
  // Gaussian elimination
  {
    // initialize sweep matrix
    std::array<double, matrix_size> sweep_matrix;
    sweep_matrix.fill(0.0);
    // prepare the linear system
    ComputeGMS(sweep_matrix,
               psi,
               cell.num_nodes,
               direction,
               cell,
               args.src_moment,
               args.groupset_start,
               idx.group_idx,
               args.num_groups,
               num_moments);
    ComputeSurfaceIntegral(sweep_matrix,
                           psi,
                           cell.num_nodes,
                           direction,
                           cell,
                           flud_data,
                           face_offset_data,
                           args.flud_stride);
    // solve
    GaussianElimination(sweep_matrix, psi, cell.num_nodes);
  }
  // post-processing
  PutPsiToFludsAndOutflow(psi,
                          direction,
                          cell,
                          flud_data,
                          face_offset_data,
                          args.flud_stride,
                          args.groupset_start + idx.group_idx,
                          cell_local_idx);
  ComputePhi(psi,
             cell.num_nodes,
             direction,
             cell,
             args.phi,
             args.groupset_start,
             idx.group_idx,
             args.num_groups,
             num_moments);
  if (psi_device)
  {
    RecordPsi(psi, cell.num_nodes, psi_device, idx, args.angleset_size, args.groupset_size);
  }
}

void
AAHSweepChunk::GPUSweep(AngleSet& angle_set)
{
  // check arguments
  if (max_num_cell_dofs_ > max_dof)
  {
    log.Log0Warning() << "GPU sweep expects max number of DOF to be " << max_dof
                      << ", but got a mesh with cell(s) requiring upto " << max_num_cell_dofs_
                      << " DOFs.\n";
    throw std::runtime_error("Max DOF exceeded. Cannot run this mesh on GPU!\n");
  }
  CALI_CXX_MARK_SCOPE("AahSweepChunk::Sweep");
  SweepKernelArgs args;
  // get mesh and quadrature data
  MeshCarrier* mesh = reinterpret_cast<MeshCarrier*>(problem_.GetCarrier(2));
  args.mesh_data = mesh->GetDevicePtr();
  QuadratureCarrier* quadrature = reinterpret_cast<QuadratureCarrier*>(groupset_.quad_carrier);
  args.quad_data = quadrature->GetDevicePtr();
  // copy source moment and destination phi data to GPU
  MemoryPinner<double>* src = reinterpret_cast<MemoryPinner<double>*>(problem_.GetPinner(0));
  src->CopyToDevice();
  args.src_moment = src->GetDevicePtr();
  MemoryPinner<double>* phi = reinterpret_cast<MemoryPinner<double>*>(problem_.GetPinner(1));
  phi->CopyToDevice();
  args.phi = phi->GetDevicePtr();
  // copy angleset data to GPU
  MemoryPinner<std::uint32_t>* directions =
    reinterpret_cast<MemoryPinner<std::uint32_t>*>(angle_set.GetMemoryPin());
  args.directions = directions->GetDevicePtr();
  std::size_t angleset_size = angle_set.GetNumAngles();
  args.angleset_size = angleset_size;
  // copy groupset data to GPU
  std::size_t groupset_size = groupset_.groups.size();
  int groupset_start = groupset_.groups.front().id;
  std::size_t num_groups = problem_.GetGroups().size();
  args.groupset_start = groupset_start;
  args.groupset_size = groupset_size;
  args.num_groups = num_groups;
  // prepare flatten FLUDs
  AAHD_FLUDS flatten_fluds(problem_, groupset_, angle_set, IsSurfaceSourceActive());
  args.flud_data = flatten_fluds.GetDevicePtr();
  args.flud_stride = flatten_fluds.GetStride();
  // prepare memory for level
  Storage<std::uint32_t>* level_storage = reinterpret_cast<Storage<std::uint32_t>*>(level_vector_);
  args.level = level_storage->GetDevicePtr();
  // allocate data for angular if required
  Storage<double> psi_storage;
  double* psi_device = nullptr;
  if (save_angular_flux_)
  {
    psi_storage = Storage<double>(max_level_size_ * angleset_size * groupset_size * max_dof);
    psi_device = psi_storage.GetDevicePtr();
  }
  // loop ever each SPDS level
  const SPDS& spds = angle_set.GetSPDS();
  const std::vector<std::vector<int>>& levelized_spls = spds.GetLevelizedLocalSubgrid();
  for (std::uint32_t level_idx = 0; const std::vector<int>& level : levelized_spls)
  {
    // copy cell local index to GPU
    level_storage->Copy(level.begin(), level.end());
    args.level_size = level.size();
    // update batch size
    std::uint32_t batch_size = level.size() * angleset_size * groupset_size;
    args.batch_size = batch_size;
    // copy FLUD of the current level to device
    flatten_fluds.CopyToDevice(level_idx);
    // launch the sweep kernel
    std::uint32_t threads_per_block = 128; // optimal value based on hardware
    std::uint32_t num_blocks = (batch_size + threads_per_block - 1) / threads_per_block;
    SweepKernel<<<num_blocks, threads_per_block>>>(args, psi_device);
    // copy FLUDs back from GPU
    flatten_fluds.CopyFromDevice(level_idx);
    // save angular flux if needed
    if (save_angular_flux_)
    {
      psi_storage.CopyFromDevice();
      crb::HostVector<double>& psi_host = psi_storage.GetHostVector();
      // loop for each cell in level, angle in angleset and group in groupset
      std::uint32_t batch_idx = 0;
      for (const int& cell_local_id : level)
      {
        Cell& cell = grid_->local_cells[cell_local_id];
        const CellMapping& cell_mapping = discretization_.GetCellMapping(cell);
        std::uint32_t cell_num_nodes = cell_mapping.GetNumNodes();
        double* cell_psi_data =
          &destination_psi_[discretization_.MapDOFLocal(cell, 0, groupset_.psi_uk_man_, 0, 0)];
        for (std::uint32_t as_ss_idx = 0; as_ss_idx < angleset_size; ++as_ss_idx)
        {
          auto direction_num = angle_set.GetAngleIndices()[as_ss_idx];
          std::size_t direction_imap = direction_num * groupset_group_stride_;
          for (std::uint32_t g_idx = 0; g_idx < groupset_size; ++g_idx)
          {
            std::size_t angle_group_imap = direction_imap + g_idx;
            for (std::uint32_t i = 0; i < cell_num_nodes; ++i)
            {
              cell_psi_data[i * groupset_angle_group_stride_ + angle_group_imap] =
                psi_host[batch_idx * max_dof + i];
            }
            ++batch_idx;
          }
        }
      }
    }
    // advance level counter
    ++level_idx;
  }
  // copy phi back to CPU
  phi->CopyFromDevice();
  // retrieve outflow
  OutflowCarrier* outflow = reinterpret_cast<OutflowCarrier*>(problem_.GetCarrier(1));
  outflow->AccumulateBack(cell_transport_views_);
  outflow->Reset();
}

void
AAHSweepChunk::CreateDeviceLevelVector()
{
  Storage<std::uint32_t>* level_vector = new Storage<std::uint32_t>(max_level_size_);
  level_vector_ = level_vector;
}

void
AAHSweepChunk::DestroyDeviceLevelVector()
{
  if (level_vector_)
  {
    Storage<std::uint32_t>* level_vector = reinterpret_cast<Storage<std::uint32_t>*>(level_vector_);
    delete level_vector;
    level_vector_ = nullptr;
  }
}

} // namespace opensn
