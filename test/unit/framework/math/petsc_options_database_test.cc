// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/petsc_utils/petsc_utils.h"
#include <gtest/gtest.h>
#include <optional>
#include <sstream>
#include <string>

using namespace opensn;

namespace
{

/// Returns whether the global database holds the option.
bool
IsSet(const std::string& name)
{
  PetscBool set = PETSC_FALSE;
  OpenSnPETScCall(PetscOptionsHasName(nullptr, nullptr, name.c_str(), &set));
  return set == PETSC_TRUE;
}

/// Returns the value of the option in the global database, or nothing for an option without one.
std::optional<std::string>
Value(const std::string& name)
{
  const char* value = nullptr;
  PetscBool set = PETSC_FALSE;
  OpenSnPETScCall(PetscOptionsFindPair(nullptr, nullptr, name.c_str(), &value, &set));
  EXPECT_TRUE(set) << name;
  if (value == nullptr)
    return std::nullopt;
  return std::string(value);
}

/// Returns whether the global database holds an option of exactly this name, without the fallback
/// of PETSc lookups.
bool
IsInDatabase(const std::string& name)
{
  char* all_options = nullptr;
  OpenSnPETScCall(PetscOptionsGetAll(nullptr, &all_options));
  std::istringstream tokens(all_options);
  OpenSnPETScCall(PetscFree(all_options));
  std::string token;
  while (tokens >> token)
    if (token == name)
      return true;
  return false;
}

} // namespace

TEST(PETScOptionsDatabaseTest, ScopeSetsAndRestoresGlobalOptions)
{
  OpenSnPETScCall(PetscOptionsSetValue(nullptr, "-pod_a_ksp_type", "cg"));
  {
    PETScOptionsDatabase database;
    database.Insert("-pod_a_ksp_type gmres -pod_a_flag -pod_a_view \"ascii:solver report.txt\"");
    database.Insert("-pod_a_ksp_type bcgs");

    // The options are not in the global database before a scope exists.
    EXPECT_EQ(Value("-pod_a_ksp_type"), "cg");
    EXPECT_FALSE(IsSet("-pod_a_flag"));
    EXPECT_FALSE(IsSet("-pod_a_view"));
    {
      const PETScOptionsDatabase::Scope scope(database);
      // The later insertion takes precedence over the earlier one and over the global value.
      EXPECT_EQ(Value("-pod_a_ksp_type"), "bcgs");
      EXPECT_TRUE(IsSet("-pod_a_flag"));
      EXPECT_EQ(Value("-pod_a_flag"), std::nullopt);
      // A quoted value keeps its space.
      EXPECT_EQ(Value("-pod_a_view"), "ascii:solver report.txt");

      // A query of the default database with a prefix, as direct solver packages make, finds the
      // option.
      PetscBool flag = PETSC_FALSE;
      PetscBool set = PETSC_FALSE;
      OpenSnPETScCall(PetscOptionsGetBool(nullptr, "pod_a_", "-flag", &flag, &set));
      EXPECT_TRUE(set);
      EXPECT_TRUE(flag);
    }
    EXPECT_EQ(Value("-pod_a_ksp_type"), "cg");
    EXPECT_FALSE(IsSet("-pod_a_flag"));
    EXPECT_FALSE(IsSet("-pod_a_view"));

    // A second scope applies the same options again.
    {
      const PETScOptionsDatabase::Scope scope(database);
      EXPECT_EQ(Value("-pod_a_ksp_type"), "bcgs");
    }
    EXPECT_EQ(Value("-pod_a_ksp_type"), "cg");
  }
  OpenSnPETScCall(PetscOptionsClearValue(nullptr, "-pod_a_ksp_type"));
}

TEST(PETScOptionsDatabaseTest, NumberedOptionIsRestoredExactly)
{
  // PETSc finds -pod_d_1_ksp_type through -pod_d_ksp_type when the former is absent. The scope must
  // not leave the former behind with the value of the latter.
  OpenSnPETScCall(PetscOptionsSetValue(nullptr, "-pod_d_ksp_type", "cg"));
  {
    PETScOptionsDatabase database;
    database.Insert("-pod_d_1_ksp_type bcgs");
    {
      const PETScOptionsDatabase::Scope scope(database);
      EXPECT_TRUE(IsInDatabase("-pod_d_1_ksp_type"));
      EXPECT_EQ(Value("-pod_d_1_ksp_type"), "bcgs");
    }
  }
  EXPECT_FALSE(IsInDatabase("-pod_d_1_ksp_type"));
  EXPECT_EQ(Value("-pod_d_1_ksp_type"), "cg");
  OpenSnPETScCall(PetscOptionsClearValue(nullptr, "-pod_d_ksp_type"));
}

TEST(PETScOptionsDatabaseTest, NestedScopesRestoreInReverseOrder)
{
  PETScOptionsDatabase outer;
  outer.Insert("-pod_b_value 1 -pod_b_outer");
  PETScOptionsDatabase inner;
  inner.Insert("-pod_b_value 2 -pod_b_inner");
  {
    const PETScOptionsDatabase::Scope outer_scope(outer);
    EXPECT_EQ(Value("-pod_b_value"), "1");
    {
      const PETScOptionsDatabase::Scope inner_scope(inner);
      EXPECT_EQ(Value("-pod_b_value"), "2");
      EXPECT_TRUE(IsSet("-pod_b_outer"));
      EXPECT_TRUE(IsSet("-pod_b_inner"));
    }
    EXPECT_EQ(Value("-pod_b_value"), "1");
    EXPECT_TRUE(IsSet("-pod_b_outer"));
    EXPECT_FALSE(IsSet("-pod_b_inner"));
  }
  EXPECT_FALSE(IsSet("-pod_b_value"));
  EXPECT_FALSE(IsSet("-pod_b_outer"));
}

TEST(PETScOptionsDatabaseTest, BlockJacobiSubSolverReadsOptions)
{
  // A diagonally dominant tridiagonal matrix.
  const PetscInt n = 4;
  Mat A = nullptr;
  OpenSnPETScCall(MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 3, nullptr, &A));
  for (PetscInt i = 0; i < n; ++i)
  {
    OpenSnPETScCall(MatSetValue(A, i, i, 4.0, INSERT_VALUES));
    if (i > 0)
      OpenSnPETScCall(MatSetValue(A, i, i - 1, -1.0, INSERT_VALUES));
    if (i < n - 1)
      OpenSnPETScCall(MatSetValue(A, i, i + 1, -1.0, INSERT_VALUES));
  }
  OpenSnPETScCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
  OpenSnPETScCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

  KSP ksp = nullptr;
  OpenSnPETScCall(KSPCreate(PETSC_COMM_SELF, &ksp));
  OpenSnPETScCall(KSPSetOperators(ksp, A, A));
  OpenSnPETScCall(KSPSetOptionsPrefix(ksp, "pod_c_"));

  PETScOptionsDatabase database;
  database.Insert("-pod_c_pc_type bjacobi -pod_c_sub_pc_type lu");
  {
    // The block solver is created during setup, so the options must be present then.
    const PETScOptionsDatabase::Scope scope(database);
    OpenSnPETScCall(KSPSetFromOptions(ksp));
    OpenSnPETScCall(KSPSetUp(ksp));
  }
  EXPECT_FALSE(IsSet("-pod_c_sub_pc_type"));

  PC pc = nullptr;
  OpenSnPETScCall(KSPGetPC(ksp, &pc));
  PCType pc_type = nullptr;
  OpenSnPETScCall(PCGetType(pc, &pc_type));
  EXPECT_STREQ(pc_type, PCBJACOBI);

  PetscInt num_blocks = 0;
  KSP* sub_ksps = nullptr;
  OpenSnPETScCall(PCBJacobiGetSubKSP(pc, &num_blocks, nullptr, &sub_ksps));
  ASSERT_EQ(num_blocks, 1);
  PC sub_pc = nullptr;
  OpenSnPETScCall(KSPGetPC(sub_ksps[0], &sub_pc));
  PCType sub_pc_type = nullptr;
  OpenSnPETScCall(PCGetType(sub_pc, &sub_pc_type));
  EXPECT_STREQ(sub_pc_type, PCLU);

  OpenSnPETScCall(KSPDestroy(&ksp));
  OpenSnPETScCall(MatDestroy(&A));
}
