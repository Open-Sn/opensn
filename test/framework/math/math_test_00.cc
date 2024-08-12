#include "framework/math/dense_vector.h"
#include "framework/math/dense_matrix.h"
#include "framework/math/sparse_matrix/sparse_matrix.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include "lua/framework/console/console.h"

using namespace opensn;

namespace unit_tests
{

ParameterBlock math_Test00(const InputParameters& params);

RegisterWrapperFunctionInNamespace(unit_tests, math_Test00, nullptr, math_Test00);

ParameterBlock
math_Test00(const InputParameters& params)
{
  opensn::log.Log() << "GOLD_BEGIN";
  // Dynamic Vector
  {
    opensn::log.Log() << "Testing DenseVector\n";

    DenseVector<double> vec(5);
    vec.Set(1.0);

    opensn::log.Log() << vec.PrintStr();
  }
  // Dynamic Matrix
  {
    opensn::log.Log() << "Testing DenseMatrix\n";
    DenseMatrix<double> mat(5, 7);
    mat.Set(1.0);

    opensn::log.Log() << mat.PrintStr();
  }

  // SparseMatrix
  {
    SparseMatrix matrix(4, 4);
    auto& mat = matrix;
    mat.Insert(0, 0, 1.0);
    mat.Insert(0, 1, 1.1);
    mat.Insert(0, 2, 1.2);
    mat.Insert(0, 3, 1.3);
    mat.Insert(1, 0, 1.9);
    mat.Insert(1, 1, 2.0);
    mat.Insert(1, 2, 2.1);
    mat.Insert(2, 1, 2.9);
    mat.Insert(2, 2, 3.0);
    mat.Insert(2, 3, 3.1);
    mat.Insert(3, 2, 3.9);
    mat.Insert(3, 3, 4.0);

    {
      auto& m = mat;
      opensn::log.Log() << "----- SparseMatrix::PrintS() -----"
                        << "\n"
                        << m.PrintStr() << "\n";

      opensn::log.Log() << "----- for (const auto& entry : m.Row(2)) -----";
      for (const auto& entry : m.Row(2))
        opensn::log.Log() << entry.row_index << " " << entry.column_index << " " << entry.value;

      opensn::log.Log() << "----- after value*2 -----";
      for (const auto& [row_index, column_index, value] : m.Row(2))
        value *= 2;

      for (const auto& entry : m.Row(2))
        opensn::log.Log() << entry.row_index << " " << entry.column_index << " " << entry.value;
    }

    opensn::log.Log() << "----- for (auto entry : matrix) -----";
    for (const auto& entry : matrix)
      opensn::log.Log() << entry.row_index << " " << entry.column_index << " " << entry.value;

    matrix.Compress();
    opensn::log.Log() << "----- after compress -----";
    for (const auto& entry : matrix)
      opensn::log.Log() << entry.row_index << " " << entry.column_index << " " << entry.value;
  }

  opensn::log.Log() << "GOLD_END";
  return ParameterBlock();
}

} //  namespace unit_tests
