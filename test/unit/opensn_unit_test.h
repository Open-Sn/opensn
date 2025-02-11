#pragma once

#include <gtest/gtest.h>

class OpenSnUnitTest : public ::testing::Test
{
public:
  void SetUp() override;
  void TearDown() override;
};
