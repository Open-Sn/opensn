// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <type_traits>
#include <cassert>
#include <array>
#include <memory>
#include <initializer_list>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <vector>
#include <string>

namespace opensn
{

template <typename T>
class NDArray
{
private:
  static constexpr size_t max_rank = 10;

public:
  /**
   * Creates an empty array.
   * \throw std::bad_alloc if memory allocation fails.
   *
   * This constructor creates an empty array and initializes the reference
   * count to one.
   */
  NDArray() noexcept : rank_(0), size_(0), storage_(nullptr), dimensions_{}, strides_{} {}

  /**
   * Creates an array with the specified number of elements in each dimension,
   * from a vector-list.
   *
   * \param dims `std::vector` list of the number of elements in each
   *             dimension.
   * \throw std::bad_alloc if memory allocation fails.
   *
   * This constructor creates an array with the specified size.
   */
  template <typename D>
  explicit NDArray(const std::vector<D>& dims)
  {
    SetDimensions(dims);
    Initialize();
  }

  /**
   * Creates an array with the specified number of elements in each dimension,
   * from a vector-list and initializes the array.
   *
   * \param dims `std::vector` list of the number of elements in each
   *             dimension.
   * \param value Initial element value.
   * \throw std::bad_alloc if memory allocation fails.
   *
   * This constructor creates an array with the specified size.
   */
  template <typename D>
  explicit NDArray(const std::vector<D>& dims, T value)
  {
    SetDimensions(dims);
    ValueInitialize(value);
  }

  /**
   * Creates an array with the specified number of elements in each dimension,
   * from an array.
   *
   * \param dims `std::array` list of the number of elements in each
   *             dimension.
   * \throw std::bad_alloc if memory allocation fails.
   *
   * This constructor creates an array with the specified size.
   */
  template <typename D, size_t N>
  explicit NDArray(const std::array<D, N>& dims)
  {
    SetDimensions(dims);
    Initialize();
  }

  /**
   * Creates an array with the specified number of elements in each dimension,
   * from an array and initializes the array.
   *
   * \param dims `std::array` list of the number of elements in each
   *             dimension.
   * \param value Initial element value.
   * \throw std::bad_alloc if memory allocation fails.
   *
   * This constructor creates an array with the specified size.
   */
  template <typename D, size_t N>
  explicit NDArray(const std::array<D, N>& dims, T value)
  {
    SetDimensions(dims);
    ValueInitialize(value);
  }

  /**
   * Creates an array with the specified number of elements in each dimension,
   * from an initializer-list.
   *
   * \param dims `std::vector` list of the number of elements in each
   *             dimension.
   * \throw std::bad_alloc if memory allocation fails.
   *
   * This constructor creates an array with the specified size.
   */
  template <typename D>
  NDArray(const std::initializer_list<D>& dims)
  {
    SetDimensions(dims);
    Initialize();
  }

  /**
   * Creates an array with the specified number of elements in each dimension,
   * from an initializer-list and initializes the array.
   *
   * \param dims `std::vector` list of the number of elements in each
   *             dimension.
   * \param value Initial element value.
   * \throw std::bad_alloc if memory allocation fails.
   *
   * This constructor creates an array with the specified size.
   */
  template <typename D>
  NDArray(const std::initializer_list<D>& dims, T value)
  {
    SetDimensions(dims);
    ValueInitialize(value);
  }

  /// Copy constructor
  NDArray(const NDArray<T>& other)
    : rank_(other.rank_),
      size_(other.size_),
      storage_(std::make_unique<T[]>(other.size_)),
      dimensions_(other.dimensions_),
      strides_(other.strides_)
  {
    std::copy(other.storage_.get(), other.storage_.get() + size_, storage_.get());
  }

  /// Move constructor
  NDArray(NDArray<T>&& other) noexcept
    : rank_(other.rank_),
      size_(other.size_),
      storage_(std::move(other.storage_)),
      dimensions_(std::move(other.dimensions_)),
      strides_(std::move(other.strides_))
  {
    other.rank_ = 0;
    other.size_ = 0;
    other.base_ = nullptr;
  }

  /// Copy assignment operator
  NDArray<T>& operator=(const NDArray<T>& other)
  {
    if (this != &other)
    {
      if (size_ != other.size_)
        storage_ = std::make_unique<T[]>(other.size_);
      std::copy(other.storage_.get(), other.storage_.get() + other.size_, storage_.get());
      rank_ = other.rank_;
      size_ = other.size_;
      dimensions_ = other.dimensions_;
      strides_ = other.strides_;
    }
    return *this;
  }

  /// Move assignment operator
  NDArray<T>& operator=(NDArray<T>&& other) noexcept
  {
    if (this != &other)
    {
      rank_ = other.rank_;
      size_ = other.size_;
      dimensions_ = std::move(other.dimensions_);
      strides_ = std::move(other.strides_);
      storage_ = std::move(other.storage_);
      other.rank_ = 0;
      other.size_ = 0;
    }
    return *this;
  }

  /**
   * Resizes the array with a vector.
   *
   * \param dims std::vector of the number of elements in each
   *             dimension.
   * \throw std::bad_alloc if memory allocation fails.
   *
   * This method resizes the array to the specified number of elements. If the
   * current size is equal to the new size, no memory allocation occurs.
   */
  template <typename D>
  void resize(const std::vector<D>& dims)
  {
    SetDimensions(dims);
    if (size_ != std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>()))
      Initialize();
  }

  /**
   * Resizes the array with an array.
   *
   * \param dims std::array of the number of elements in each
   *             dimension.
   * \throw std::bad_alloc if memory allocation fails.
   *
   * This method resizes the array to the specified number of elements. If the
   * current size is equal to the new size, no memory allocation occurs.
   */
  template <typename D, size_t N>
  void resize(const std::array<D, N>& dims)
  {
    SetDimensions(dims);
    if (size_ != std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>()))
      Initialize();
  }

  /**
   * Resizes the array with an initializer_list.
   *
   * \param dims std::initializer_list of the number of elements in each
   *             dimension.
   * \throw std::bad_alloc if memory allocation fails.
   *
   * This method resizes the array to the specified number of elements. If the
   * current size is equal to the new size, no memory allocation occurs.
   */
  template <typename D>
  void resize(const std::initializer_list<D>& dims)
  {
    SetDimensions(dims);
    if (size_ != std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>()))
      Initialize();
  }

  /**
   * Accesses the specified element for an array with N dimensions.
   *
   * \param args Indices for each dimension.
   * \return Read/write reference to the element.
   */
  template <typename... Args>
  inline __attribute__((always_inline)) constexpr T& operator()(Args... args) noexcept
  {
    return storage_[ComputeIndex(args...)];
  }

  /**
   * Accesses the specified element for an array with N dimensions with bounds checking.
   *
   * \param args Indices for each dimension.
   * \return Read/write reference to the element.
   */
  template <typename... Args>
  constexpr T& at(Args... args)
  {
    size_t indices[]{static_cast<size_t>(args)...};
    for (size_t i = 0; i < rank_; ++i)
    {
      if (indices[i] >= dimensions_[i])
        throw std::out_of_range("Index out of bounds.");
    }
    return storage_[ComputeIndex(args...)];
  }

  /**
   * Accesses the specified element for an array with N dimensions.
   *
   * \param args Indices for each dimension.
   * \return Read reference to the element.
   */
  template <typename... Args>
  inline __attribute__((always_inline)) constexpr const T& operator()(Args... args) const noexcept
  {
    return storage_[ComputeIndex(args...)];
  }

  /**
   * Accesses the specified element for an array with N dimensions with bounds checking.
   *
   * \param args Indices for each dimension.
   * \return Read reference to the element.
   */
  template <typename... Args>
  constexpr const T& at(Args... args) const
  {
    size_t indices[]{static_cast<size_t>(args)...};
    for (size_t i = 0; i < rank_; ++i)
    {
      if (indices[i] >= dimensions_[i])
        throw std::out_of_range("Index out of bounds.");
    }
    return storage_[ComputeIndex(args...)];
  }

  /// Returns an iterator pointing to the beginning of the array.
  inline constexpr T* begin() const noexcept { return storage_.get(); }

  /// Returns a constant iterator pointing to the beginning of the array.
  inline constexpr const T* cbegin() const noexcept { return storage_.get(); }

  /// Returns an iterator pointing to the end of the array.
  inline constexpr T* end() const noexcept { return storage_.get() + size_; }

  /// Returns a constant iterator pointing to the end of the array.
  inline constexpr const T* cend() const noexcept { return storage_.get() + size_; }

  /// Returns the number of elements in the array.
  inline constexpr size_t size() const noexcept { return size_; }

  /// Returns true if the array has no elements
  inline constexpr bool empty() const noexcept { return size_ == 0; }

  /// Returns a pointer to the underlying array data.
  inline constexpr T* data() const noexcept { return storage_.get(); }

  /// Returns the rank of the array.
  inline constexpr size_t rank() const noexcept { return rank_; }

  /// Returns the dimension of the array.
  inline std::vector<size_t> dimension() const noexcept
  {
    return std::vector<size_t>(dimensions_.begin(), dimensions_.begin() + rank_);
  }

  /// Sets each element of the array to the specified value.
  inline void set(T value) noexcept { std::fill(storage_.get(), storage_.get() + size_, value); }

  /**
   * Returns a linear index to the specified element with safety checks.
   *
   * \param args The indices of the desired element.
   * \throw std::invalid_argument if the number of arguments are incorrect and
   * std::out_of_range if one of the dimension-indices are out of range.
   * \return Linear index to the specified element.
   */
  template <typename... Args>
  size_t MapNDtoLin(Args... args) const
  {
    if (sizeof...(args) != rank_)
    {
      throw std::invalid_argument("Number of arguments " + std::to_string(sizeof...(args)) +
                                  " not equal to rank " + std::to_string(rank_));
    }

    const size_t N = rank_;
    size_t indices[]{static_cast<size_t>(args)...};
    for (size_t i = 0; i < N; ++i)
    {
      if (indices[i] >= dimensions_[i])
      {
        throw std::out_of_range("Index " + std::to_string(i) + " out of range " +
                                std::to_string(indices[i]) + " must be <" +
                                std::to_string(dimensions_[i]));
      }
    }

    return ComputeIndex(args...);
  }

  /**
   * Swap the contents of this array with another array.
   *
   * \param other The array to swap with.
   */
  inline void swap(NDArray<T>& other) noexcept
  {
    std::swap(rank_, other.rank_);
    std::swap(size_, other.size_);
    std::swap(storage_, other.storage_);
    std::swap(dimensions_, other.dimensions_);
    std::swap(strides_, other.strides_);
  }

private:
  template <typename D>
  void SetDimensions(const std::vector<D>& dims)
  {
    rank_ = dims.size();
    if (rank_ > max_rank)
      throw std::invalid_argument("Rank exceeds the maximum allowed value of 10.");
    std::copy(dims.begin(), dims.end(), dimensions_.begin());
  }

  template <typename D, size_t N>
  void SetDimensions(const std::array<D, N>& dims)
  {
    rank_ = N;
    static_assert(N <= max_rank, "Rank exceeds the maximum allowed value of 10.");
    std::copy(dims.begin(), dims.end(), dimensions_.begin());
  }

  template <typename D>
  void SetDimensions(const std::initializer_list<D>& dims)
  {
    rank_ = dims.size();
    if (rank_ > max_rank)
      throw std::invalid_argument("Rank exceeds the maximum allowed value of 10.");
    std::copy(dims.begin(), dims.end(), dimensions_.begin());
  }

  void Initialize()
  {
    size_ = 1;
    strides_[rank_ - 1] = 1;
    for (size_t i = rank_; i-- > 0;)
    {
      size_ *= dimensions_[i];
      if (i > 0)
        strides_[i - 1] = strides_[i] * dimensions_[i];
    }

    storage_ = std::make_unique<T[]>(size_);
  }

  void ValueInitialize(T value)
  {
    Initialize();
    std::fill_n(storage_.get(), size_, value);
  }

  template <typename... Args>
  inline __attribute__((always_inline)) constexpr size_t ComputeIndex(Args... args) const noexcept
  {
    size_t indices[]{static_cast<size_t>(args)...};

    if constexpr (sizeof...(args) == 1)
      return indices[0];
    else if constexpr (sizeof...(args) == 2)
      return indices[0] * strides_[0] + indices[1] * strides_[1];
    else if constexpr (sizeof...(args) == 3)
      return indices[0] * strides_[0] + indices[1] * strides_[1] + indices[2] * strides_[2];
    else if constexpr (sizeof...(args) == 4)
    {
      return indices[0] * strides_[0] + indices[1] * strides_[1] + indices[2] * strides_[2] +
             indices[3] * strides_[3];
    }
    else if constexpr (sizeof...(args) == 5)
    {
      return indices[0] * strides_[0] + indices[1] * strides_[1] + indices[2] * strides_[2] +
             indices[3] * strides_[3] + indices[4] * strides_[4];
    }

    size_t index = 0;
    for (size_t i = 0; i < rank_; ++i)
      index += indices[i] * strides_[i];
    return index;
  }

private:
  size_t rank_;
  size_t size_;
  std::unique_ptr<T[]> storage_;
  std::array<size_t, max_rank> dimensions_;
  std::array<size_t, max_rank> strides_;
};

} // namespace opensn
