#include "framework/math/dynamic_vector.h"
#include "framework/math/dynamic_matrix.h"
#include "framework/math/SparseMatrix/chi_math_sparse_matrix.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include "framework/console/chi_console.h"

namespace chi_unit_tests
{

chi::ParameterBlock chi_math_Test00(const chi::InputParameters& params);

RegisterWrapperFunction(/*namespace_name=*/chi_unit_tests,
                        /*name_in_lua=*/chi_math_Test00,
                        /*syntax_function=*/nullptr,
                        /*actual_function=*/chi_math_Test00);

chi::ParameterBlock
chi_math_Test00(const chi::InputParameters& params)
{
  Chi::log.Log() << "GOLD_BEGIN";
  // Dynamic Vector
  {
    Chi::log.Log() << "Testing chi_math::DynamicVector\n";

    chi_math::DynamicVector<double> vec(5, 1.0);

    Chi::log.Log() << vec.PrintStr();
  }
  // Dynamic Matrix
  {
    Chi::log.Log() << "Testing chi_math::DynamicMatrix\n";
    chi_math::DynamicMatrix<double> mat(5, 7, 1.0);

    Chi::log.Log() << mat.PrintStr();
  }

  // SparseMatrix
  {
    chi_math::SparseMatrix matrix(4, 4);
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
      Chi::log.Log() << "----- SparseMatrix::PrintS() -----"
                     << "\n"
                     << m.PrintStr() << "\n";

      Chi::log.Log() << "----- for (const auto& entry : m.Row(2)) -----";
      for (const auto& entry : m.Row(2))
        Chi::log.Log() << entry.row_index << " " << entry.column_index << " " << entry.value;

      Chi::log.Log() << "----- after value*2 -----";
      for (const auto& [row_index, column_index, value] : m.Row(2))
        value *= 2;

      for (const auto& entry : m.Row(2))
        Chi::log.Log() << entry.row_index << " " << entry.column_index << " " << entry.value;
    }

    Chi::log.Log() << "----- for (auto entry : matrix) -----";
    for (const auto& entry : matrix)
      Chi::log.Log() << entry.row_index << " " << entry.column_index << " " << entry.value;

    matrix.Compress();
    Chi::log.Log() << "----- after compress -----";
    for (const auto& entry : matrix)
      Chi::log.Log() << entry.row_index << " " << entry.column_index << " " << entry.value;
  }

  Chi::log.Log() << "GOLD_END";
  return chi::ParameterBlock();
}

} // namespace chi_unit_tests
