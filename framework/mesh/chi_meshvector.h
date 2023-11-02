#ifndef CHI_MESH_VECTOR3_H
#define CHI_MESH_VECTOR3_H
#include <iostream>
#include <cmath>
#include <sstream>

namespace chi_mesh
{
struct TensorRank2Dim3;
}

// namespace chi_mesh
//{
//=============================================== General 3D vector structure
/**General 3 element vector structure.
 * \author Jan
 */
struct chi_mesh::Vector3
{
  double x; ///< Element-0
  double y; ///< Element-1
  double z; ///< Element-2

  /**Default constructor. Initialized as all zeros.*/
  Vector3()
  {
    x = 0.0;
    y = 0.0;
    z = 0.0;
  }

  /**Constructor where single element is initialized \f$ x[z]=a \f$.*/
  explicit Vector3(double a)
  {
    x = a;
    y = 0.0;
    z = 0.0;
  }

  /**Constructor where \f$ x=a\f$ and \f$ y=b \f$. */
  explicit Vector3(double a, double b)
  {
    x = a;
    y = b;
    z = 0.0;
  }

  /**Constructor where \f$ \vec{x}=[a,b,c] \f$.*/
  explicit Vector3(double a, double b, double c)
  {
    x = a;
    y = b;
    z = c;
  }

  /**Constructor where \f$ \vec{x}=\{a,b,c\} \f$.*/
  Vector3(std::initializer_list<double> list)
  {
    if (not empty(list))
    {
      std::vector<double> vec = list;
      for (size_t i = 0; ((i < 3) and (i < vec.size())); ++i)
      {
        if (i == 0) x = vec[i];
        if (i == 1) y = vec[i];
        if (i == 2) z = vec[i];
      }
    }
  }

  /**Constructor where \f$ \vec{x}=\{a,b,c\} \f$.*/
  explicit Vector3(const std::vector<double>& list)
  {
    if (not empty(list))
    {
      std::vector<double> vec = list;
      for (size_t i = 0; ((i < 3) and (i < vec.size())); ++i)
      {
        if (i == 0) x = vec[i];
        if (i == 1) y = vec[i];
        if (i == 2) z = vec[i];
      }
    }
  }

  /**Copy constructor.*/
  Vector3(const Vector3& that)
  {
    this->x = that.x;
    this->y = that.y;
    this->z = that.z;
  }

  /**Assignment operator.*/
  Vector3& operator=(const Vector3& that)
  {
    this->x = that.x;
    this->y = that.y;
    this->z = that.z;

    return *this;
  }

  Vector3& operator=(std::initializer_list<double> list)
  {
    if (not empty(list))
    {
      std::vector<double> vec = list;
      for (size_t i = 0; ((i < 3) and (i < vec.size())); ++i)
      {
        this->operator()(i) = vec[i];
      }
    }

    return *this;
  }

  Vector3& operator=(const std::vector<double>& list)
  {
    if (not empty(list))
    {
      std::vector<double> vec = list;
      for (size_t i = 0; ((i < 3) and (i < vec.size())); ++i)
      {
        this->operator()(i) = vec[i];
      }
    }

    return *this;
  }

  //============================================= Addition
  /**Component-wise addition of two vectors.
   * \f$ \vec{w} = \vec{x} + \vec{y} \f$*/
  Vector3 operator+(const Vector3& that) const
  {
    Vector3 newVector;
    newVector.x = this->x + that.x;
    newVector.y = this->y + that.y;
    newVector.z = this->z + that.z;

    return newVector;
  }

  /**In-place component-wise addition of two vectors.
   * \f$ \vec{x} = \vec{x} + \vec{y} \f$*/
  Vector3& operator+=(const Vector3& that)
  {
    this->x += that.x;
    this->y += that.y;
    this->z += that.z;

    return *this;
  }

  /**Component-wise shift by scalar-value.
   * \f$ \vec{w} = \vec{x} + \alpha \f$*/
  Vector3 Shifted(const double value) const
  {
    Vector3 newVector;
    newVector.x = this->x + value;
    newVector.y = this->y + value;
    newVector.z = this->z + value;

    return newVector;
  }

  /**In-place component-wise shift by scalar-value.
   * \f$ \vec{x} = \vec{x} + \alpha \f$*/
  Vector3& Shift(const double value)
  {
    this->x += value;
    this->y += value;
    this->z += value;

    return *this;
  }

  //============================================= Subtraction
  /**Component-wise subtraction.
   * \f$ \vec{w} = \vec{x} - \vec{y} \f$*/
  Vector3 operator-(const Vector3& that) const
  {
    Vector3 newVector;
    newVector.x = this->x - that.x;
    newVector.y = this->y - that.y;
    newVector.z = this->z - that.z;

    return newVector;
  }

  /**In-place component-wise subtraction.
   * \f$ \vec{x} = \vec{x} - \vec{y} \f$*/
  Vector3& operator-=(const Vector3& that)
  {
    this->x -= that.x;
    this->y -= that.y;
    this->z -= that.z;

    return *this;
  }

  //============================================= Multiplication
  /**Vector component-wise multiplication by scalar.
   * \f$ \vec{w} = \vec{x} \alpha \f$*/
  Vector3 operator*(const double value) const
  {
    Vector3 newVector;
    newVector.x = this->x * value;
    newVector.y = this->y * value;
    newVector.z = this->z * value;

    return newVector;
  }

  /**Vector in-place component-wise multiplication by scalar.
   * \f$ \vec{x} = \vec{x} \alpha \f$*/
  Vector3& operator*=(const double value)
  {
    this->x *= value;
    this->y *= value;
    this->z *= value;

    return *this;
  }

  /**Vector component-wise multiplication.
   * \f$ w_i = x_i y_i \f$*/
  Vector3 operator*(const Vector3& that) const
  {
    Vector3 newVector;
    newVector.x = this->x * that.x;
    newVector.y = this->y * that.y;
    newVector.z = this->z * that.z;

    return newVector;
  }

  /**Vector in-place component-wise multiplication.
   * \f$ x_i = x_i y_i \f$*/
  Vector3& operator*=(const Vector3& that)
  {
    this->x *= that.x;
    this->y *= that.y;
    this->z *= that.z;

    return *this;
  }

  //============================================= Division
  /**Vector component-wise division by scalar.
   * \f$ w_i = \frac{x_i}{\alpha} \f$*/
  Vector3 operator/(const double value) const
  {
    Vector3 newVector;
    newVector.x = this->x / value;
    newVector.y = this->y / value;
    newVector.z = this->z / value;

    return newVector;
  }

  /**Vector in-place component-wise division by scalar.
   * \f$ x_i = \frac{x_i}{\alpha} \f$*/
  Vector3& operator/=(const double value)
  {
    this->x /= value;
    this->y /= value;
    this->z /= value;

    return *this;
  }

  /**Vector component-wise division.
   * \f$ w_i = \frac{x_i}{y_i} \f$*/
  Vector3 operator/(const Vector3& that) const
  {
    Vector3 newVector;
    newVector.x = this->x / that.x;
    newVector.y = this->y / that.y;
    newVector.z = this->z / that.z;

    return newVector;
  }

  /**Vector in-place component-wise division.
   * \f$ x_i = \frac{x_i}{y_i} \f$*/
  Vector3& operator/=(const Vector3& that)
  {
    this->x /= that.x;
    this->y /= that.y;
    this->z /= that.z;

    return *this;
  }

  //============================================= Element access
  /**Returns a copy of the value at the given index.*/
  double operator[](const size_t i) const
  {
    if (i == 0) return this->x;
    else if (i == 1)
      return this->y;
    else if (i == 2)
      return this->z;

    return 0.0;
  }

  /**Returns a reference of the value at the given index.*/
  double& operator()(const size_t i)
  {
    if (i == 0) return this->x;
    else if (i == 1)
      return this->y;
    else if (i == 2)
      return this->z;

    return this->x;
  }

  //============================================= Tensor product
  // Defined in chi_mesh_utilities.cc
  chi_mesh::TensorRank2Dim3 OTimes(const Vector3& that) const;

  //============================================= Tensor dot product
  // Defined in chi_mesh_utilities.cc
  Vector3 Dot(const chi_mesh::TensorRank2Dim3& that) const;

  //============================================= Operations
  /**Vector cross-product.
   * \f$ \vec{w} = \vec{x} \times \vec{y} \f$*/
  Vector3 Cross(const Vector3& that) const
  {
    Vector3 newVector;
    newVector.x = this->y * that.z - this->z * that.y;
    newVector.y = this->z * that.x - this->x * that.z;
    newVector.z = this->x * that.y - this->y * that.x;

    return newVector;
  }

  /**Vector dot-product.
   * \f$ \vec{w} = \vec{x} \bullet \vec{y} \f$ */
  double Dot(const Vector3& that) const
  {
    double value = 0.0;
    value += this->x * that.x;
    value += this->y * that.y;
    value += this->z * that.z;

    return value;
  }

  /**Computes the L2-norm of the vector. Otherwise known as the length of
   * a 3D vector.*/
  double Norm() const
  {
    double value = 0.0;
    value += this->x * this->x;
    value += this->y * this->y;
    value += this->z * this->z;

    value = sqrt(value);

    return value;
  }

  /**Computes the square of the L2-norm of the vector. This eliminates the
   * usage of the square root and is therefore less expensive that a proper
   * L2-norm. Useful if only comparing distances.*/
  double NormSquare() const
  {
    double value = 0.0;
    value += this->x * this->x;
    value += this->y * this->y;
    value += this->z * this->z;

    return value;
  }

  /**Normalizes the vector in-place.
   * \f$ \vec{x} = \frac{\vec{x}}{||x||_2} \f$*/
  void Normalize()
  {
    double norm = this->Norm();

    x /= norm;
    y /= norm;
    z /= norm;
  }

  /**Returns a normalized version of the vector.
   * \f$ \vec{w} = \frac{\vec{x}}{||x||_2} \f$*/
  Vector3 Normalized() const
  {
    double norm = this->Norm();

    Vector3 newVector;
    newVector.x = this->x / norm;
    newVector.y = this->y / norm;
    newVector.z = this->z / norm;

    return newVector;
  }

  /**Returns a vector v^* where each element is inverted provided
   * that it is greater than the given tolerance, otherwise the offending entry
   * is set to 0.0.
   * \f$ w_i = \frac{1.0}{x_i} \f$*/
  Vector3 InverseZeroIfSmaller(const double tol) const
  {
    Vector3 newVector;
    newVector.x = (std::fabs(this->x) > tol) ? 1.0 / this->x : 0.0;
    newVector.y = (std::fabs(this->y) > tol) ? 1.0 / this->y : 0.0;
    newVector.z = (std::fabs(this->z) > tol) ? 1.0 / this->z : 0.0;

    return newVector;
  }

  /**Returns a vector v^* where each element is inverted provided
   * that it is greater than the given tolerance, otherwise the offending entry
   * is set to 1.0.
   * \f$ w_i = \frac{1.0}{x_i} \f$*/
  Vector3 InverseOneIfSmaller(const double tol) const
  {
    Vector3 newVector;
    newVector.x = (std::fabs(this->x) > tol) ? 1.0 / this->x : 1.0;
    newVector.y = (std::fabs(this->y) > tol) ? 1.0 / this->y : 1.0;
    newVector.z = (std::fabs(this->z) > tol) ? 1.0 / this->z : 1.0;

    return newVector;
  }

  /**Returns a vector v^* where each element is inverted without any
   * check for division by zero.
   * \f$ w_i = \frac{1.0}{x_i} \f$*/
  Vector3 Inverse() const
  {
    Vector3 newVector;
    double dx_inv = 1.0 / this->x;
    double dy_inv = 1.0 / this->y;
    double dz_inv = 1.0 / this->z;

    return newVector;
  }

  /**Prints the vector to std::cout.*/
  void Print() const
  {
    std::cout << this->x << " ";
    std::cout << this->y << " ";
    std::cout << this->z;
  }

  friend std::ostream& operator<<(std::ostream& out, Vector3& v)
  {
    out << "[" << v.x << " " << v.y << " " << v.z << "]";

    return out;
  }

  /**Deprecated. Prints the vector to a string and then returns the string.*/
  std::string PrintS() const // TODO: Deprecated
  {
    std::stringstream out;
    out << "[" << x << " " << y << " " << z << "]";

    return out.str();
  }

  /**Prints the vector to a string and then returns the string.*/
  std::string PrintStr() const
  {
    std::stringstream out;
    out << "[" << x << " " << y << " " << z << "]";

    return out.str();
  }

  static size_t Size() { return 3; }
};

// The following functions are defined in chi_mesh_utilities.cc
// Left multiplcation by scalar
chi_mesh::Vector3 operator*(double value, const chi_mesh::Vector3& that);

//}//namespace chi_mesh

#endif // CHI_MESH_VECTOR3_H
