#pragma once

#include <gtest/gtest.h>
#include "framework/mesh/mesh_continuum/mesh_continuum.h"

using namespace opensn;

class OpenSnUnitTest : public ::testing::Test
{
public:
  void SetUp() override;
  void TearDown() override;
};
