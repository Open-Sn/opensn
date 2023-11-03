#pragma once

#include <cstddef>
#include <type_traits>
#include <cassert>
#include <vector>
#include <array>
#include <stdexcept>
#include <string>

namespace chi_data_types
{

template <typename T>
class NDArray
{
private:
  size_t rank_;
  size_t* dimensions_;
  size_t* strides_;
  size_t size_;
  T* base_;

  template <bool...>
  struct bool_pack
  {
  };

  template <class... U>
  using conjunction = std::is_same<bool_pack<true, U::value...>, bool_pack<U::value..., true>>;

  template <typename... U>
  using AllIntegral = typename conjunction<std::is_integral<U>...>::type;

public:
  /** Creates an array with the specified number of elements in each dimension,
   *  from a vector-list.
   *
   *  \param dims `std::vector` list of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template <typename D>
  explicit NDArray(const std::vector<D>& dims)
    : rank_(0), dimensions_(nullptr), strides_(nullptr), size_(0), base_(nullptr)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");
    const size_t N = dims.size();

    rank_ = N;

    strides_ = new size_t[N];
    dimensions_ = new size_t[N];

    // Populate dimensions
    for (size_t i = 0; i < N; ++i)
      dimensions_[i] = dims[i];

    // Populate strides and size
    size_ = 1;
    for (size_t i = 0; i < N; ++i)
    {
      assert(dimensions_[i] > 0);
      size_ *= dimensions_[i];
      strides_[i] = 1;
      for (size_t j = i + 1; j < N; ++j)
        strides_[i] *= dimensions_[j];
    }

    base_ = new T[size_];
    for (size_t i = 0; i < size_; ++i)
      base_[i] = 0.0;
  }

  /** Creates an array with the specified number of elements in each dimension,
   *  from an array.
   *
   *  \param dims `std::array` list of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template <typename D, size_t N>
  explicit NDArray(const std::array<D, N>& dims)
    : rank_(0), dimensions_(nullptr), strides_(nullptr), size_(0), base_(nullptr)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have array of integral types.");
    rank_ = N;

    strides_ = new size_t[N];
    dimensions_ = new size_t[N];

    // Populate dimensions
    for (size_t i = 0; i < N; ++i)
      dimensions_[i] = dims[i];

    // Populate strides and size
    size_ = 1;
    for (size_t i = 0; i < N; ++i)
    {
      assert(dimensions_[i] > 0);
      size_ *= dimensions_[i];
      strides_[i] = 1;
      for (size_t j = i + 1; j < N; ++j)
        strides_[i] *= dimensions_[j];
    }

    base_ = new T[size_];
    for (size_t i = 0; i < size_; ++i)
      base_[i] = 0.0;
  }

  /** Creates an array with the specified number of elements in each dimension,
   *  from an initializer-list.
   *
   *  \param dims `std::vector` list of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template <typename D>
  NDArray(const std::initializer_list<D>& dims)
    : rank_(0), dimensions_(nullptr), strides_(nullptr), size_(0), base_(nullptr)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");
    const size_t N = dims.size();

    rank_ = N;

    strides_ = new size_t[N];
    dimensions_ = new size_t[N];

    // Populate dimensions
    {
      size_t i = 0;
      for (size_t val : dims)
      {
        dimensions_[i] = val;
        ++i;
        if (i >= N) break;
      }
    }

    // Populate strides and size
    size_ = 1;
    for (size_t i = 0; i < N; ++i)
    {
      assert(dimensions_[i] > 0);
      size_ *= dimensions_[i];
      strides_[i] = 1;
      for (size_t j = i + 1; j < N; ++j)
        strides_[i] *= dimensions_[j];
    }

    base_ = new T[size_];
    for (size_t i = 0; i < size_; ++i)
      base_[i] = 0.0;
  }

  /** Creates an array with the specified number of elements in each dimension,
   *  from a vector. Each entry in the array is assigned the designated value.
   *
   *  \param dims `std::vector` list of the number of elements in each
   *              dimension.
   *  \param value The value to assing to each element.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template <typename D>
  explicit NDArray(const std::vector<D>& dims, T value)
    : rank_(0), dimensions_(nullptr), strides_(nullptr), size_(0), base_(nullptr)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");
    const size_t N = dims.size();

    rank_ = N;

    strides_ = new size_t[N];
    dimensions_ = new size_t[N];

    // Populate dimensions
    for (size_t i = 0; i < N; ++i)
      dimensions_[i] = dims[i];

    // Populate strides and size
    size_ = 1;
    for (size_t i = 0; i < N; ++i)
    {
      assert(dimensions_[i] > 0);
      size_ *= dimensions_[i];
      strides_[i] = 1;
      for (size_t j = i + 1; j < N; ++j)
        strides_[i] *= dimensions_[j];
    }

    base_ = new T[size_];
    for (size_t i = 0; i < size_; ++i)
      base_[i] = value;
  }

  /** Creates an array with the specified number of elements in each dimension,
   *  from an array. Each entry in the array is assigned the designated value.
   *
   *  \param dims `std::array` list of the number of elements in each
   *              dimension.
   *  \param value The value to assing to each element.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template <typename D, size_t N>
  explicit NDArray(const std::array<D, N>& dims, T value)
    : rank_(0), dimensions_(nullptr), strides_(nullptr), size_(0), base_(nullptr)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");
    rank_ = N;

    strides_ = new size_t[N];
    dimensions_ = new size_t[N];

    // Populate dimensions
    for (size_t i = 0; i < N; ++i)
      dimensions_[i] = dims[i];

    // Populate strides and size
    size_ = 1;
    for (size_t i = 0; i < N; ++i)
    {
      assert(dimensions_[i] > 0);
      size_ *= dimensions_[i];
      strides_[i] = 1;
      for (size_t j = i + 1; j < N; ++j)
        strides_[i] *= dimensions_[j];
    }

    base_ = new T[size_];
    for (size_t i = 0; i < size_; ++i)
      base_[i] = value;
  }

  /** Creates an array with the specified number of elements in each dimension,
   *  from an initializer-list. Each entry in the array is assigned the
   *  designated value.
   *
   *  \param dims `std::initializer` list of the number of elements in each
   *              dimension.
   *  \param value The value to assing to each element.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an array with the specified size.
   */
  template <typename D>
  NDArray(const std::initializer_list<D>& dims, T value)
    : rank_(0), dimensions_(nullptr), strides_(nullptr), size_(0), base_(nullptr)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");
    const size_t N = dims.size();

    rank_ = N;

    strides_ = new size_t[N];
    dimensions_ = new size_t[N];

    // Populate dimensions
    {
      size_t i = 0;
      for (size_t val : dims)
      {
        dimensions_[i] = val;
        ++i;
        if (i >= N) break;
      }
    }

    // Populate strides and size
    size_ = 1;
    for (size_t i = 0; i < N; ++i)
    {
      assert(dimensions_[i] > 0);
      size_ *= dimensions_[i];
      strides_[i] = 1;
      for (size_t j = i + 1; j < N; ++j)
        strides_[i] *= dimensions_[j];
    }

    base_ = new T[size_];
    for (size_t i = 0; i < size_; ++i)
      base_[i] = value;
  }

  /** Creates an empty array.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This constructor creates an empty array and initializes the reference
   * count to one.
   */
  NDArray() : rank_(0), dimensions_(nullptr), strides_(nullptr), size_(0), base_(nullptr) {}

  /** Copy construct from another array.
   *  \param other The array to copy.
   */
  NDArray(NDArray<T> const& other)
    : rank_(other.rank_),
      dimensions_(nullptr),
      strides_(nullptr),
      size_(other.size_),
      base_(nullptr)
  {

    const size_t N = rank_;
    strides_ = new size_t[N];
    dimensions_ = new size_t[N];
    base_ = new T[size_];

    for (size_t i = 0; i < N; ++i)
    {
      dimensions_[i] = other.dimensions_[i];
      strides_[i] = other.strides_[i];
    }

    for (size_t i = 0; i < size_; ++i)
      base_[i] = other.base_[i];
  }

  /** Assign from another array.
   *  \param other The array to copy.
   */
  NDArray<T>& operator=(NDArray<T> const& other)
  {
    NDArray<T>(other).swap(*this);
    return *this;
  }

  /**Move constructor.*/
  NDArray(NDArray<T>&& other) noexcept
    : rank_(std::move(other.rank_)),
      dimensions_(std::move(other.dimensions_)),
      strides_(std::move(other.strides_)),
      size_(std::move(other.size_)),
      base_(std::move(other.base_))
  {
  }

  /**Deleted move assignment operator*/
  NDArray<T>& operator=(NDArray<T>&&) = delete;

  /**Sets a value to all the items in the array.*/
  void set(T value)
  {
    for (size_t i = 0; i < size_; ++i)
      base_[i] = value;
  }

  /** Swap the contents of this array with another array.
   *  \param other The array to swap with.
   */
  void swap(NDArray<T>& other)
  {
    std::swap(rank_, other.rank_);
    std::swap(dimensions_, other.dimensions_);
    std::swap(strides_, other.strides_);
    std::swap(size_, other.size_);
    std::swap(base_, other.base_);
  }

  /** Returns the number of elements in the array.*/
  size_t size() const noexcept { return size_; }

  /** Returns true if the array has no elements.*/
  bool empty() const noexcept { return size_ == 0; }

  /** Returns an iterator pointing to the beginning of the array.*/
  T* begin() const noexcept { return base_; }

  /** Returns a constant iterator pointing to the beginning of the array.*/
  const T* cbegin() const noexcept { return base_; }

  /** Returns an iterator pointing to the end of the array.*/
  T* end() const noexcept { return base_ + size_; }

  /** Returns a constant iterator pointing to the end of the array*/
  const T* cend() const noexcept { return base_ + size_; }

  /** Returns a pointer to the underlying array data.*/
  T* data() const noexcept { return base_; }

  /** Returns the rank of the array.*/
  size_t rank() const noexcept { return rank_; }

  /**Returns the dimension of the array.*/
  std::vector<size_t> dimension() const
  {
    std::vector<size_t> dim(rank_, 0);
    for (size_t i = 0; i < rank_; ++i)
      dim[i] = dimensions_[i];

    return dim;
  }

  /** Resizes the array with a vector.
   *
   *  \param dims std::vector of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This method resizes the array to the specified number of elements. If the
   *  current size is equal to the new size, no memory allocation occurs.
   */
  template <typename D>
  void resize(const std::vector<D>& dims)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");

    delete[] dimensions_;
    dimensions_ = nullptr;
    delete[] strides_;
    strides_ = nullptr;
    delete[] base_;
    base_ = nullptr;

    NDArray<T>(dims).swap(*this);
  }

  /** Resizes the array with an array.
   *
   *  \param dims std::array of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This method resizes the array to the specified number of elements. If the
   *  current size is equal to the new size, no memory allocation occurs.
   */
  template <typename D, size_t N>
  void resize(const std::array<D, N>& dims)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");

    delete[] dimensions_;
    dimensions_ = nullptr;
    delete[] strides_;
    strides_ = nullptr;
    delete[] base_;
    base_ = nullptr;

    NDArray<T>(dims).swap(*this);
  }

  /** Resizes the array with an initializer_list.
   *
   *  \param dims std::initializer_list of the number of elements in each
   *              dimension.
   *  \throw std::bad_alloc if memory allocation fails.
   *
   *  This method resizes the array to the specified number of elements. If the
   *  current size is equal to the new size, no memory allocation occurs.
   */
  template <typename D>
  void resize(const std::initializer_list<D>& dims)
  {
    static_assert(std::is_integral<D>::value,
                  "NDArray dims argument must have vector of integral types.");

    delete[] dimensions_;
    dimensions_ = nullptr;
    delete[] strides_;
    strides_ = nullptr;
    delete[] base_;
    base_ = nullptr;

    NDArray<T>(dims).swap(*this);
  }

  /** Accesses the specified element.
   *  \param args The indices of the desired element.
   *  \return Read/write reference to the element.
   */
  template <typename... Args>
  T& operator()(Args... args) noexcept
  {
    static_assert(AllIntegral<Args...>::value,
                  "NDArray::operator[]: All parameters must be of integral type");

    size_t indices[]{static_cast<size_t>(args)...};

    const size_t N = rank_;

    T* address = base_ + indices[N - 1];
    for (size_t i = 0; i < N - 1; ++i)
      address += strides_[i] * indices[i];
    return *(address);
  }

  /** Accesses the specified element.
   *  \param args The indices of the desired element.
   *  \return Read reference to the element.
   */
  template <typename... Args>
  T const& operator()(Args... args) const noexcept
  {
    static_assert(AllIntegral<Args...>::value,
                  "NDArray::operator[]: All parameters must be of integral type");

    size_t indices[]{static_cast<size_t>(args)...};

    const size_t N = rank_;

    T* address = base_ + indices[N - 1];
    for (size_t i = 0; i < N - 1; ++i)
      address += strides_[i] * indices[i];
    return *(address);
  }

  /** Accesses the specified element with safety checks.
   *
   *  \param args The indices of the desired element.
   *  \throw std::invalid_argument if the number of arguments are incorrect and
   *  std::out_of_range if one of the dimension-indices are out of range.
   *  \return Read/write reference to the element.
   */
  template <typename... Args>
  T& at(Args... args) noexcept
  {
    static_assert(AllIntegral<Args...>::value,
                  "NDArray::at(): All parameters must be of integral type");

    if (sizeof...(args) != rank_)
      throw std::invalid_argument("NDArray::at(): Number of arguments " +
                                  std::to_string(sizeof...(args)) + " not equal to rank " +
                                  std::to_string(rank_));

    const size_t N = rank_;
    size_t indices[]{static_cast<size_t>(args)...};
    for (size_t i = 0; i < N; ++i)
      if (indices[i] >= dimensions_[i])
        throw std::out_of_range("NDArray::at(): Index " + std::to_string(i) + " out of range " +
                                std::to_string(indices[i]) + " must be <" +
                                std::to_string(dimensions_[i]));

    T* address = base_ + indices[N - 1];
    for (size_t i = 0; i < N - 1; ++i)
      address += strides_[i] * indices[i];
    return *(address);
  }

  /** Returns a linear index to the specified element with safety checks.
   *
   *  \param args The indices of the desired element.
   *  \throw std::invalid_argument if the number of arguments are incorrect and
   *  std::out_of_range if one of the dimension-indices are out of range.
   *  \return Linear index to the specified element.
   */
  template <typename... Args>
  size_t MapNDtoLin(Args... args) const
  {
    static_assert(AllIntegral<Args...>::value,
                  "NDArray::at(): All parameters must be of integral type");

    if (sizeof...(args) != rank_)
      throw std::invalid_argument("NDArray::at(): Number of arguments " +
                                  std::to_string(sizeof...(args)) + " not equal to rank " +
                                  std::to_string(rank_));

    const size_t N = rank_;
    size_t indices[]{static_cast<size_t>(args)...};
    for (size_t i = 0; i < N; ++i)
      if (indices[i] >= dimensions_[i])
        throw std::out_of_range("NDArray::at(): Index " + std::to_string(i) + " out of range " +
                                std::to_string(indices[i]) + " must be <" +
                                std::to_string(dimensions_[i]));

    size_t index = indices[N - 1];
    for (size_t i = 0; i < N - 1; ++i)
      index += strides_[i] * indices[i];
    return index;
  }

  /** Deletes the array.
   *
   *  The destructor deletes the underlying array data.
   */
  ~NDArray()
  {
    delete[] dimensions_;
    dimensions_ = nullptr;
    delete[] strides_;
    strides_ = nullptr;
    delete[] base_;
    base_ = nullptr;
  }
};

} // namespace chi_data_types
