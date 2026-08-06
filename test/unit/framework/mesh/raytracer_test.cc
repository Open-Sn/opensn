// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "gtest/gtest.h"
#include "test/unit/common/mesh_builders.h"
#include "framework/mesh/raytrace/raytracer.h"
#include "framework/runtime.h"
#include <cmath>

using namespace opensn;

namespace
{

const Cell&
FindCellByCentroid(const std::shared_ptr<MeshContinuum>& grid, const Vector3& centroid)
{
  for (const auto& cell : grid->GetLocalCells())
    if (cell->centroid.AbsoluteEquals(centroid, 1.0e-12))
      return *cell;

  throw std::logic_error("Failed to find test cell by centroid.");
}

} // namespace

TEST(RayTracerTest, TraceRayHandlesFaceAndVertexStarts)
{
  if (opensn::mpi_comm.size() != 1)
    GTEST_SKIP() << "Ray tracer endpoint test is serial only.";

  const auto grid = BuildSquareMesh(2.0, 2, -1.0);
  RayTracer ray_tracer(grid);

  const auto& lower_right_cell = FindCellByCentroid(grid, Vector3(0.5, -0.5, 0.0));
  const Vector3 face_start_initial(0.0, -0.5, 0.0);
  Vector3 face_start = face_start_initial;
  Vector3 face_direction(1.0, 0.0, 0.0);
  auto face_trace = ray_tracer.TraceRay(lower_right_cell, face_start, face_direction);

  EXPECT_FALSE(face_trace.particle_lost);
  EXPECT_GT(face_trace.distance_to_surface, 0.0);
  EXPECT_NEAR(face_trace.distance_to_surface, 1.0, 1.0e-10);
  EXPECT_TRUE(face_start.AbsoluteEquals(face_start_initial, 1.0e-14));
  EXPECT_TRUE(face_trace.pos_f.AbsoluteEquals(Vector3(1.0, -0.5, 0.0), 1.0e-10));
  EXPECT_NEAR((face_trace.pos_f - face_start_initial).Norm(), 1.0, 1.0e-10);

  const auto& upper_right_cell = FindCellByCentroid(grid, Vector3(0.5, 0.5, 0.0));
  const Vector3 vertex_start_initial(0.0, 0.0, 0.0);
  Vector3 vertex_start = vertex_start_initial;
  Vector3 vertex_direction = Vector3(1.0, 1.0, 0.0).Normalized();
  auto vertex_trace = ray_tracer.TraceRay(upper_right_cell, vertex_start, vertex_direction);

  EXPECT_FALSE(vertex_trace.particle_lost);
  EXPECT_GT(vertex_trace.distance_to_surface, 0.0);
  EXPECT_NEAR(vertex_trace.distance_to_surface, std::sqrt(2.0), 1.0e-10);
  EXPECT_TRUE(vertex_start.AbsoluteEquals(vertex_start_initial, 1.0e-14));
  EXPECT_TRUE(vertex_trace.pos_f.AbsoluteEquals(Vector3(1.0, 1.0, 0.0), 1.0e-10));
  EXPECT_NEAR((vertex_trace.pos_f - vertex_start_initial).Norm(), std::sqrt(2.0), 1.0e-10);
}

TEST(RayTracerTest, TraceRayHandlesEdgeStart3D)
{
  if (opensn::mpi_comm.size() != 1)
    GTEST_SKIP() << "Ray tracer endpoint test is serial only.";

  const auto grid = BuildBoxMesh(2.0, 2, -1.0);
  RayTracer ray_tracer(grid);

  const auto& upper_right_front_cell = FindCellByCentroid(grid, Vector3(0.5, 0.5, 0.5));
  const Vector3 edge_start_initial(0.0, 0.0, 0.5);
  Vector3 edge_start = edge_start_initial;
  Vector3 edge_direction = Vector3(1.0, 1.0, 0.0).Normalized();
  auto edge_trace = ray_tracer.TraceRay(upper_right_front_cell, edge_start, edge_direction);

  EXPECT_FALSE(edge_trace.particle_lost);
  EXPECT_GT(edge_trace.distance_to_surface, 0.0);
  EXPECT_NEAR(edge_trace.distance_to_surface, std::sqrt(2.0), 1.0e-10);
  EXPECT_TRUE(edge_start.AbsoluteEquals(edge_start_initial, 1.0e-14));
  EXPECT_TRUE(edge_trace.pos_f.AbsoluteEquals(Vector3(1.0, 1.0, 0.5), 1.0e-10));
  EXPECT_NEAR((edge_trace.pos_f - edge_start_initial).Norm(), std::sqrt(2.0), 1.0e-10);
}
