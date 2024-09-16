#include "lua/framework/console/console.h"
#include "framework/parameters/input_parameters.h"
#include "framework/math/dense_matrix.h"
#include "framework/logging/log.h"
#include "framework/runtime.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock
DenseMatrixTest(const InputParameters&)
{
  OpenSnLogicalErrorIf(opensn::mpi_comm.size() != 1, "Requires 1 processor");

  opensn::log.Log() << "GOLD_BEGIN";
  DenseMatrix<double> a(2, 4);
  opensn::log.LogAll() << "a:size = " << a.Rows() << ", " << a.Columns() << std::endl;

  a(0, 0) = 1;
  a(0, 1) = 2;
  a(0, 2) = 3;
  a(0, 4) = 4;
  a(1, 0) = -1;
  a(1, 1) = -2;
  a(1, 2) = -3;
  a(1, 4) = -4;
  opensn::log.LogAll() << "a:values = " << std::endl;
  opensn::log.LogAll() << a.PrintStr() << std::endl;

  Vector<double> r0(4);
  r0(0) = 10;
  r0(1) = 11;
  r0(2) = 12;
  r0(3) = 13;
  a.SetRow(0, r0);
  opensn::log.LogAll() << "a:values = " << std::endl;
  opensn::log.LogAll() << a.PrintStr() << std::endl;

  SwapRows(a, 0, 1);
  opensn::log.LogAll() << "a:values = " << std::endl;
  opensn::log.LogAll() << a.PrintStr() << std::endl;

  DenseMatrix<double> diag(3, 3, 0.);
  diag.SetDiagonal(12.);
  opensn::log.LogAll() << "diag:values = " << std::endl;
  opensn::log.LogAll() << diag.PrintStr() << std::endl;

  auto a_trans = Transpose(a);
  opensn::log.LogAll() << "a_trans:values = " << std::endl;
  opensn::log.LogAll() << a_trans.PrintStr() << std::endl;

  auto a_mult = Mult(a, 2.);
  opensn::log.LogAll() << "a_mult:values = " << std::endl;
  opensn::log.LogAll() << a_mult.PrintStr() << std::endl;

  Vector<double> v0(4);
  v0(0) = 1;
  v0(1) = -1;
  v0(2) = 3;
  v0(3) = -2;
  auto a_mult_vec = Mult(a, v0);
  opensn::log.LogAll() << "a_mult_vec:values = " << std::endl;
  opensn::log.LogAll() << a_mult_vec.PrintStr() << std::endl;

  DenseMatrix<double> b(4, 3);
  b(0, 0) = -1;
  b(0, 1) = 0;
  b(0, 2) = 2;
  b(1, 0) = -4;
  b(1, 1) = 0;
  b(1, 2) = -2;
  b(2, 0) = 3;
  b(2, 1) = 5;
  b(2, 2) = 8;
  b(3, 0) = 1;
  b(3, 1) = -1;
  b(3, 2) = 0;
  auto ab = Mult(a, b);
  opensn::log.LogAll() << "ab:values = " << std::endl;
  opensn::log.LogAll() << ab.PrintStr() << std::endl;

  auto apb = Add(a, a_mult);
  opensn::log.LogAll() << "a+b:values = " << std::endl;
  opensn::log.LogAll() << apb.PrintStr() << std::endl;

  auto amb = Subtract(a, a_mult);
  opensn::log.LogAll() << "a-b:values = " << std::endl;
  opensn::log.LogAll() << amb.PrintStr() << std::endl;

  auto a_scaled = a;
  Scale(a_scaled, 0.5);
  opensn::log.LogAll() << "a_scaled:values = " << std::endl;
  opensn::log.LogAll() << a_scaled.PrintStr() << std::endl;

  a_scaled = Scaled(a, 0.5);
  opensn::log.LogAll() << "a_scaled:values = " << std::endl;
  opensn::log.LogAll() << a_scaled.PrintStr() << std::endl;

  auto d0 = Determinant(diag);
  opensn::log.LogAll() << "d0 = " << d0 << std::endl;

  opensn::log.Log() << "GOLD_END";

  return ParameterBlock();
}

RegisterWrapperFunctionInNamespace(unit_tests, DenseMatrixTest, nullptr, DenseMatrixTest);

} // namespace unit_tests
