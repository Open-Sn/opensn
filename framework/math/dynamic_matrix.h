#pragma once

#include <vector>
#include <stdexcept>
#include <sstream>

#include "framework/math/math.h"
#include "framework/math/dynamic_vector.h"

namespace opensn
{
template <class NumberFormat>
class DynamicMatrix;

/**General dynamic matrix utility.*/
template <class NumberFormat>
class DynamicMatrix
{
  typedef std::pair<size_t, size_t> MatDim;

public:
  std::vector<std::vector<NumberFormat>> elements_;

  /**Default constructor. Does nothing.*/
  DynamicMatrix()
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DynamicMatrix.");
  }

  /**Constructor with number of entries. Value defaults.*/
  DynamicMatrix(size_t Nrows, size_t Ncols) : elements_(Nrows, std::vector<NumberFormat>(Ncols))
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DynamicMatrix.");
  }

  /**Constructor with number of entries and default value.*/
  DynamicMatrix(size_t Nrows, size_t Ncols, NumberFormat value)
    : elements_(Nrows, std::vector<NumberFormat>(Ncols, value))
  {
    static_assert(std::is_floating_point<NumberFormat>::value,
                  "Only floating point number formats are "
                  "supported for DynamicMatrix.");
  }

  /**Copy constructor.*/
  DynamicMatrix(const DynamicMatrix& other) { elements_ = other.elements_; }

  /**Assignment operator.*/
  DynamicMatrix& operator=(const DynamicMatrix& other)
  {
    elements_ = other.elements_;
    return *this;
  }

  /**Move constructor.*/
  DynamicMatrix(DynamicMatrix&& other) noexcept { elements_ = std::move(other.elements_); }

  /**Move assignment operator.*/
  DynamicMatrix& operator=(DynamicMatrix&& other) noexcept
  {
    elements_ = std::move(other.elements_);
    return *this;
  }

  /**Constructor with vector.*/
  explicit DynamicMatrix(const std::vector<std::vector<double>>& in) { elements_ = in; }

  /**Copy constructor with vector.*/
  DynamicMatrix& operator=(const std::vector<std::vector<double>>& in)
  {
    elements_ = in;
    return *this;
  }

  /**Constructor with vector.*/
  DynamicMatrix(std::initializer_list<std::initializer_list<NumberFormat>> in)
  {
    elements_.clear();
    for (auto& v : in)
      elements_.push_back(v);
  }

  /**Copy constructor with vector.*/
  DynamicMatrix& operator=(std::initializer_list<std::initializer_list<NumberFormat>> in)
  {
    elements_.clear();
    for (auto& v : in)
      elements_.push_back(v);
    return *this;
  }

  // Element access
  std::vector<NumberFormat>& operator[](size_t i) { return elements_[i]; }

  std::vector<NumberFormat>& at(size_t i) { return elements_.at(i); }

  std::vector<NumberFormat>& back() { return elements_.back(); }

  std::vector<NumberFormat>& front() { return elements_.front(); }

  std::vector<NumberFormat>* data() { return elements_.data(); }

  void clear() { elements_.clear(); }

  void resize(size_t Nrows, size_t Ncols)
  {
    elements_.resize(Nrows);
    for (auto& row : elements_)
      row.resize(Ncols);
  }

  void resize(size_t Nrows, size_t Ncols, const NumberFormat& val)
  {
    elements_.resize(Nrows);
    for (auto& row : elements_)
    {
      row.resize(Ncols, val);
      for (auto& entry : row)
        entry = val;
    }
  }

  void reserve(size_t Nrows) { elements_.reserve(Nrows); }

  void push_back(const std::vector<NumberFormat>& val) { elements_.push_back(val); }
  void pop_back() { elements_.pop_back(); }

  bool empty() const noexcept { return elements_.empty(); }

  // Iterator access
  typename std::vector<std::vector<NumberFormat>>::iterator begin() { return elements_.begin(); }

  typename std::vector<std::vector<NumberFormat>>::iterator end() { return elements_.end(); }

  size_t size() const { return elements_.size(); }

  MatDim Dimensions() const
  {
    if (elements_.empty())
      return {0, 0};
    else
      return {elements_.size(), elements_[0].size()};
  }

private:
  static void bounds_check_rows_cols(const MatDim a, const MatDim b)
  {
    if ((a.first != b.first) or (a.second != b.second))
      throw std::length_error("Mismatched square sizes of DynamicMatrix");
  }

  static void bounds_check_colsA_rowsB(const MatDim a, const MatDim b)
  {
    if (a.first != b.second)
      throw std::length_error("Mismatched matrix A rows with matrix B cols"
                              " in DynamicMatrix");
  }

public:
  /**Component-wise addition of two matrices.
   * \f$ \vec{w} = \vec{x} + \vec{y} \f$*/
  DynamicMatrix operator+(const DynamicMatrix& rhs) const
  {
    auto dim = Dimensions();
    bounds_check_rows_cols(dim, rhs.Dimensions());
    DynamicMatrix<NumberFormat> newVector(dim.first, dim.second, 0.0);
    for (int i = 0; i < dim.first; ++i)
      for (int j = 0; j < dim.second; ++j)
        newVector.elements_[i][j] = elements_[i][j] + rhs.elements_[i][j];

    return newVector;
  }

  /**In-place component-wise addition of two vectors.
   * \f$ \vec{x} = \vec{x} + \vec{y} \f$*/
  DynamicMatrix& operator+=(const DynamicMatrix& rhs)
  {
    auto dim = Dimensions();
    bounds_check_rows_cols(dim, rhs.Dimensions());
    for (int i = 0; i < dim.first; ++i)
      for (int j = 0; j < dim.second; ++j)
        elements_[i][j] = elements_[i][j] + rhs.elements_[i][j];

    return *this;
  }

  /**Component-wise subtraction.
   * \f$ \vec{w} = \vec{x} - \vec{y} \f$*/
  DynamicMatrix operator-(const DynamicMatrix& rhs) const
  {
    auto dim = Dimensions();
    bounds_check_rows_cols(dim, rhs.Dimensions());
    DynamicMatrix<NumberFormat> newVector(dim.first, dim.second, 0.0);
    for (int i = 0; i < dim.first; ++i)
      for (int j = 0; j < dim.second; ++j)
        newVector.elements_[i][j] = elements_[i][j] - rhs.elements_[i][j];

    return newVector;
  }

  /**In-place component-wise subtraction.
   * \f$ \vec{x} = \vec{x} - \vec{y} \f$*/
  DynamicMatrix& operator-=(const DynamicMatrix& rhs)
  {
    auto dim = Dimensions();
    bounds_check_rows_cols(dim, rhs.Dimensions());
    for (int i = 0; i < dim.first; ++i)
      for (int j = 0; j < dim.second; ++j)
        elements_[i][j] = elements_[i][j] - rhs.elements_[i][j];

    return *this;
  }

  /**Vector component-wise multiplication by scalar.
   * \f$ \vec{w} = \vec{x} \alpha \f$*/
  DynamicMatrix operator*(const NumberFormat value) const
  {
    auto dim = Dimensions();
    DynamicMatrix<NumberFormat> newVector(dim.first, dim.second);
    for (int i = 0; i < dim.first; ++i)
      for (int j = 0; j < dim.second; ++j)
        newVector.elements_[i][j] = elements_[i][j] * value;

    return newVector;
  }

  /**Vector in-place component-wise multiplication by scalar.
   * \f$ \vec{x} = \vec{x} \alpha \f$*/
  DynamicMatrix& operator*=(const NumberFormat value)
  {
    auto dim = Dimensions();
    for (int i = 0; i < dim.first; ++i)
      for (int j = 0; j < dim.second; ++j)
        elements_[i][j] *= value;

    return *this;
  }

  /** Matrix-Matrix multiplication */
  DynamicMatrix operator*(const DynamicMatrix& rhs)
  {
    auto dimA = Dimensions();
    auto dimB = rhs.Dimensions();
    bounds_check_colsA_rowsB(dimA, dimB);

    DynamicMatrix<NumberFormat> newMatrix(dimA.first, dimB.second);
    unsigned int istar = 0;
    unsigned int jstar = 0;
    for (unsigned int i = 0; i < dimA.first; ++i)
    {
      for (unsigned int j = 0; j < dimA.first; ++j)
      {
        NumberFormat value = 0.0;
        for (unsigned int k = 0; k < dimA.second; ++k)
          value += elements_[i][k] * rhs.elements_[k][j];

        newMatrix[istar][jstar] = value;
        ++jstar;
      }
      ++istar;
      jstar = 0;
    } // for i

    return newMatrix;
  }

  /** Matrix-Vector multiplication */
  DynamicVector<NumberFormat> operator*(const DynamicVector<NumberFormat>& V)
  {
    auto dimA = Dimensions();
    auto dimV = V.size();

    if (dimA.second != dimV)
      throw std::length_error("Mismatched matrix vector sizes in"
                              " matrix-vector multiplication: DynamicMatrix");

    DynamicVector<NumberFormat> newV(dimA.first);
    unsigned int k = 0;
    for (unsigned int i = 0; i < dimA.first; ++i)
    {
      NumberFormat value = 0.0;
      for (unsigned int j = 0; j < dimA.second; ++j)
        value += elements_[i][j] * V[j];
      newV[k] = value;
      ++k;
    }

    return newV;
  }

  /**Vector component-wise division by scalar.
   * \f$ w_i = \frac{x_i}{\alpha} \f$*/
  DynamicMatrix operator/(const NumberFormat value) const
  {
    auto dim = Dimensions();
    DynamicMatrix<NumberFormat> newVector(dim.first, dim.second);
    for (int i = 0; i < dim.first; ++i)
      for (int j = 0; j < dim.second; ++j)
        newVector.elements_[i][j] = elements_[i][j] / value;

    return newVector;
  }

  /**Vector in-place component-wise division by scalar.
   * \f$ x_i = \frac{x_i}{\alpha} \f$*/
  DynamicMatrix& operator/=(const NumberFormat value)
  {
    auto dim = Dimensions();
    for (int i = 0; i < dim.first; ++i)
      for (int j = 0; j < dim.second; ++j)
        elements_[i][j] /= value;

    return *this;
  }

  /**Obtains the inverse with Gauss-Elimination.*/
  DynamicMatrix Inverse() const
  {
    auto inv_elems = InverseGEPivoting(elements_);

    return DynamicMatrix<NumberFormat>(inv_elems);
  }

  /** Set the diagonal using a vector.*/
  void SetDiagonal(DynamicVector<NumberFormat>& V)
  {
    auto dimA = Dimensions();
    auto dimV = V.size();

    if ((dimA.first != dimV) or (dimA.second != dimV))
    {
      throw std::length_error("Mismatched matrix vector sizes in"
                              " matrix-vector diagonal assignment: DynamicMatrix");
    }

    for (int i = 0; i < dimA.first; ++i)
      elements_[i][i] = V[i];
  }

  /** Set the diagonal using value.*/
  void SetDiagonal(NumberFormat val)
  {
    auto dimA = Dimensions();

    for (int i = 0; i < dimA.first; ++i)
      elements_[i][i] = val;
  }

  /**Prints the matrix to a string and then returns the string.*/
  std::string PrintStr() const
  {
    auto dim = Dimensions();
    std::stringstream out;

    for (int i = 0; i < dim.first; ++i)
    {
      for (int j = 0; j < (dim.second - 1); ++j)
        out << elements_[i][j] << " ";
      out << elements_[i][dim.second - 1];

      if (i < (dim.first - 1))
        out << "\n";
    }

    return out.str();
  }
};

/**Multiplication by a scalar from the left.*/
template <class NumberFormat>
DynamicMatrix<NumberFormat>
operator*(const double value, DynamicMatrix<NumberFormat>& that)
{
  auto dim = that.Dimensions();
  DynamicMatrix<NumberFormat> newMatrix(dim.first, dim.second);
  for (int i = 0; i < dim.first; ++i)
    for (int j = 0; j < dim.second; ++j)
      newMatrix.elements_[i][j] = that.elements_[i][j] * value;

  return newMatrix;
}

} // namespace opensn
