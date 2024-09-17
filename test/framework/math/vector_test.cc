#include "lua/framework/console/console.h"
#include "framework/parameters/input_parameters.h"
#include "framework/math/vector.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock
VectorTest(const InputParameters&)
{
  OpenSnLogicalErrorIf(opensn::mpi_comm.size() != 1, "Requires 1 processor");

  opensn::log.Log() << "GOLD_BEGIN";

  Vector<double> a(3);
  opensn::log.LogAll() << "a:size = " << a.Rows() << std::endl;
  a.Resize(4);
  opensn::log.LogAll() << "a:size = " << a.Rows() << std::endl;

  a.Set(123.);
  opensn::log.LogAll() << "a:values = " << a.PrintStr() << std::endl;

  a.Scale(0.1);
  opensn::log.LogAll() << "a:values = " << a.PrintStr() << std::endl;

  Vector<double> b(2);
  b(0) = 3;
  b(1) = 4;
  opensn::log.LogAll() << "b:values = " << b.PrintStr() << std::endl;
  auto b_l2_norm = Vec2Norm(b);

  Vector<double> c(2);
  c(0) = 5;
  c(1) = 2;
  auto d = Add(b, c);
  opensn::log.LogAll() << "b+c:values = " << d.PrintStr() << std::endl;

  auto e = Subtract(b, c);
  opensn::log.LogAll() << "b-c:values = " << e.PrintStr() << std::endl;

  auto dp0 = Dot(b, c);
  opensn::log.LogAll() << "b dot c:values = " << dp0 << std::endl;

  auto b_scaled = Scaled(b, 4.);
  opensn::log.LogAll() << "b_scaled:values = " << b_scaled.PrintStr() << std::endl;
  Scale(b_scaled, 2.);
  opensn::log.LogAll() << "b_scaled:values = " << b_scaled.PrintStr() << std::endl;

  //

  auto b_normalized = b.Normalized();
  opensn::log.LogAll() << "b_normalized = " << b_normalized.PrintStr() << std::endl;

  Vector<double> added(2);
  added.Set(1.);
  added.Add(b);
  opensn::log.LogAll() << "added = " << added.PrintStr() << std::endl;

  Vector<double> subtracted(2);
  subtracted.Set(2.);
  subtracted.Subtract(b);
  opensn::log.LogAll() << "subtracted = " << subtracted.PrintStr() << std::endl;

  Vector<double> scaled(2);
  scaled(0) = 5;
  scaled(1) = 2;
  scaled.Scale(2);
  opensn::log.LogAll() << "scaled = " << scaled.PrintStr() << std::endl;

  auto scaled2 = scaled.Scaled(0.5);
  opensn::log.LogAll() << "scaled2 = " << scaled2.PrintStr() << std::endl;

  auto dp1 = b.Dot(c);
  opensn::log.LogAll() << "b dot c = " << dp1 << std::endl;

  opensn::log.Log() << "GOLD_END";

  return ParameterBlock();
}

RegisterWrapperFunctionInNamespace(unit_tests, VectorTest, nullptr, VectorTest);

} // namespace unit_tests
