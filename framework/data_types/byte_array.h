#pragma once

#include <cstddef>
#include <vector>
#include <stdexcept>
#include <string>

namespace chi_data_types
{
class ByteArray
{
protected:
  std::vector<std::byte> raw_data_;
  size_t offset_ = 0;

public:
  ByteArray() = default;
  explicit ByteArray(const size_t raw_data_size) : raw_data_(raw_data_size) {}
  explicit ByteArray(std::vector<std::byte>&& raw_data) : raw_data_(std::move(raw_data)) {}
  explicit ByteArray(const std::vector<std::byte>& raw_data) : raw_data_(raw_data) {}

  /**Uses the template type T to convert an associated value (of type T)
   * to a sub-array of std::bytes and adds it to the internal byte-array.
   *
   * The template type T must support sizeof.*/
  template <typename T>
  void Write(const T& value)
  {
    const size_t num_bytes = sizeof(T);
    const std::byte* value_byte_array = reinterpret_cast<const std::byte*>(&value);

    raw_data_.insert(raw_data_.end(), value_byte_array, value_byte_array + num_bytes);
  }

  /**Uses the template type `T` to convert `sizeof(T)` number of bytes to
   * a value of type `T`. The bytes are pulled from the internal byte-array
   * starting at the internal address `m_offset` which can be viewed with a
   * call to `Offset()`. This offset is incremented by the amount of bytes
   * used and can be reset with a call to `Seek(0)`.
   *
   * Bounds-checking is performed by checking if the internal byte array
   * has the required number of bytes available. If this check fails then
   * this call will return an `out_of_range` exception.*/
  template <typename T>
  T Read()
  {
    const size_t num_bytes = sizeof(T);
    if ((offset_ + num_bytes - 1) >= raw_data_.size())
      throw std::out_of_range(
        std::string("ByteArray reading error. ") + " Typename: " + std::string(typeid(T).name()) +
        " m_offset: " + std::to_string(offset_) + " size: " + std::to_string(raw_data_.size()) +
        " num_bytes to read: " + std::to_string(num_bytes));

    T value = *reinterpret_cast<T*>(&raw_data_[offset_]);
    offset_ += num_bytes;

    return value;
  }

  /**Uses the template type T to convert sizeof(T) number of bytes to
   * a value of type T. The bytes are pulled from the internal byte-array
   * starting at the internal address specified by the argument "address".
   * An optional argument next_address can be used to return the location
   * of the value after the value read.
   * This offset is the given "address" argument incremented by the amount
   * of bytes used.
   *
   * Bounds-checking is performed by checking if the internal byte array
   * has the required number of bytes available. If this check fails then
   * this call will return an `out_of_range` exception.*/
  template <typename T>
  T Read(const size_t address, size_t* next_address = nullptr) const
  {
    const size_t num_bytes = sizeof(T);
    if ((address + num_bytes - 1) >= raw_data_.size())
      throw std::logic_error(
        std::string("ByteArray reading error. ") + " Typename: " + std::string(typeid(T).name()) +
        " address: " + std::to_string(address) + " size: " + std::to_string(raw_data_.size()) +
        " num_bytes to read: " + std::to_string(num_bytes));

    T value = *reinterpret_cast<const T*>(&raw_data_[address]);
    if (next_address != nullptr) *next_address = address + num_bytes;

    return value;
  }

  /**Appends a `ByteArray` to the current internal byte array.*/
  void Append(const ByteArray& other_raw)
  {
    auto& master = raw_data_;

    const auto& slave = other_raw.Data();
    master.insert(master.end(), slave.begin(), slave.end());
  }

  /**Appends bytes from a `std::vector<std::byte>` to the internal
   * byte array.*/
  void Append(const std::vector<std::byte>& other_raw)
  {
    auto& master = raw_data_;

    const auto& slave = other_raw;
    master.insert(master.end(), slave.begin(), slave.end());
  }

  /**Clears the internal byte array and resets the address offset to zero.*/
  void Clear()
  {
    raw_data_.clear();
    offset_ = 0;
  }

  /**Moves the address marker to the supplied address.*/
  void Seek(const size_t address = 0) { offset_ = address; }
  /**Returns the internal address marker position.*/
  size_t Offset() const { return offset_; }
  /**Determines if the internal address marker is beyond the internal
   * byte array.*/
  bool EndOfBuffer() const { return offset_ >= raw_data_.size(); }
  /**Returns the current size of the internal byte array.*/
  size_t Size() const { return raw_data_.size(); }

  /**Returns a reference to the internal byte array.*/
  std::vector<std::byte>& Data() { return raw_data_; }
  /**Returns a const reference of the internal byte array.*/
  const std::vector<std::byte>& Data() const { return raw_data_; }
};
} // namespace chi_data_types


