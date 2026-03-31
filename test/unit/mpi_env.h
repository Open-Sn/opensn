#pragma once

#include <gtest/gtest.h>

class MPIEnvironment : public ::testing::Environment
{
public:
  MPIEnvironment() = default;

  void SetUp() override;
  void TearDown() override;
};
