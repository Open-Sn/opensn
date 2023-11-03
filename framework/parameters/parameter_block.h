#pragma once

#include "data_types/varying.h"

#include <memory>
#include <vector>
#include <string>

namespace chi
{

enum class ParameterBlockType
{
  INVALID_VALUE = 0,
  BOOLEAN = 1,
  FLOAT = 3,
  STRING = 4,
  INTEGER = 5,
  ARRAY = 98,
  BLOCK = 99
};

std::string ParameterBlockTypeName(ParameterBlockType type);

class ParameterBlock;

/**A ParameterBlock is a conceptually simple data structure that supports
 * a hierarchy of primitive parameters. There really are just 4 member variables
 * on a ParameterBlock object, they are 1) the type (as an enum), 2) the
 * name of the block, 3) a pointer to a value (which can only be a primitive
 * type), and 4) a vector of child parameters.
 *
 * If a ParameterBlock has a primitive type, i.e., BOOLEAN, FLOAT, STRING, or
 * INTEGER, then the value_ptr will contain a pointer to the value of a
 * primitive type. Otherwise, for types ARRAY and BLOCK, the ParameterBlock
 * will not have a value_ptr and instead the vector member will contain
 * sub-parameters.*/
class ParameterBlock
{
private:
  ParameterBlockType type_ = ParameterBlockType::BLOCK;
  std::string name_;
  std::unique_ptr<chi_data_types::Varying> value_ptr_ = nullptr;
  std::vector<ParameterBlock> parameters_;
  std::string error_origin_scope_ = "Unknown Scope";

public:
  /**Sets the name of the block.*/
  void SetBlockName(const std::string& name);

public:
  // Helpers
  template <typename T>
  struct IsBool
  {
    static constexpr bool value = std::is_same_v<T, bool>;
  };
  template <typename T>
  struct IsFloat
  {
    static constexpr bool value = std::is_floating_point_v<T>;
  };
  template <typename T>
  struct IsString
  {
    static constexpr bool value = std::is_same_v<T, std::string> or std::is_same_v<T, const char*>;
  };
  template <typename T>
  struct IsInteger
  {
    static constexpr bool value = std::is_integral_v<T> and not std::is_same_v<T, bool>;
  };

  // Constructors
  /**Constructs an empty parameter block with the given name and type BLOCK.*/
  explicit ParameterBlock(const std::string& name = "");

  /**Derived type constructor*/
  template <typename T>
  ParameterBlock(const std::string& name, const std::vector<T>& array)
    : type_(ParameterBlockType::ARRAY), name_(name)
  {
    size_t k = 0;
    for (const T& value : array)
      AddParameter(std::to_string(k++), value);
  }

  /**Constructs one of the fundamental types.*/
  template <typename T>
  explicit ParameterBlock(const std::string& name, T value) : name_(name)
  {
    constexpr bool is_supported =
      IsBool<T>::value or IsFloat<T>::value or IsString<T>::value or IsInteger<T>::value;

    static_assert(is_supported, "Value type not supported for parameter block");

    if (IsBool<T>::value) type_ = ParameterBlockType::BOOLEAN;
    if (IsFloat<T>::value) type_ = ParameterBlockType::FLOAT;
    if (IsString<T>::value) type_ = ParameterBlockType::STRING;
    if (IsInteger<T>::value) type_ = ParameterBlockType::INTEGER;

    value_ptr_ = std::make_unique<chi_data_types::Varying>(value);
  }

  /**Copy constructor*/
  ParameterBlock(const ParameterBlock& other);

  /**Copy assignment operator*/
  ParameterBlock& operator=(const ParameterBlock& other);

  /**Move constructor*/
  ParameterBlock(ParameterBlock&& other) noexcept;

  /**Move assignment operator*/
  ParameterBlock& operator=(ParameterBlock&& other) noexcept;

  // Accessors
  ParameterBlockType Type() const;
  /**Returns true if the parameter block comprises a single value of any of
   * the types BOOLEAN, FLOAT, STRING, INTEGER.*/
  bool IsScalar() const;
  /**Returns a string version of the type.*/
  std::string TypeName() const;
  std::string Name() const;
  const chi_data_types::Varying& Value() const;
  size_t NumParameters() const;
  /**Returns the sub-parameters of this block.*/
  const std::vector<ParameterBlock>& Parameters() const;
  /**Returns whether or not the block has a value. If this block has
   * sub-parameters it should not have a value. This is a good way to
   * check if the block is actually a single value.*/
  bool HasValue() const;

  // Mutators
  /**Changes the block type to array, making it accessible via integer
   * keys.*/
  void ChangeToArray();

  /**Sets a string to be displayed alongside exceptions that give some
   * notion of the origin of the error.*/
  void SetErrorOriginScope(const std::string& scope);

  /**Gets a string that allows error messages to print the scope of an
   * error.*/
  std::string GetErrorOriginScope() const { return error_origin_scope_; }

  // Requirements
  /**Checks that the block is of the given type. If it is not it
   * will throw an exception `std::logic_error`.*/
  void RequireBlockTypeIs(ParameterBlockType type) const;
  void RequireParameterBlockTypeIs(const std::string& param_name, ParameterBlockType type) const
  {
    GetParam(param_name).RequireBlockTypeIs(type);
  }
  /**Check that the parameter with the given name exists otherwise
   * throws a `std::logic_error`.*/
  void RequireParameter(const std::string& param_name) const;

  // utilities
  /**Adds a parameter to the sub-parameter list.*/
  void AddParameter(ParameterBlock block);
  /**Makes a ParameterBlock and adds it to the sub-parameters list.*/
  template <typename T>
  void AddParameter(const std::string& name, T value)
  {
    AddParameter(ParameterBlock(name, value));
  }

  /**Sorts the sub-parameter list according to name. This is useful
   * for regression testing.*/
  void SortParameters();

  /**Returns true if a parameter with the specified name is in the
   * list of sub-parameters. Otherwise, false.*/
  bool Has(const std::string& param_name) const;

  /**Gets a parameter by name.*/
  ParameterBlock& GetParam(const std::string& param_name);
  /**Gets a parameter by index.*/
  ParameterBlock& GetParam(size_t index);

  /**Gets a parameter by name.*/
  const ParameterBlock& GetParam(const std::string& param_name) const;
  /**Gets a parameter by index.*/
  const ParameterBlock& GetParam(size_t index) const;

public:
  /**Returns the value of the parameter.*/
  template <typename T>
  T GetValue() const
  {
    if (value_ptr_ == nullptr)
      throw std::logic_error(error_origin_scope_ + std::string(__PRETTY_FUNCTION__) +
                             ": Value not available for block type " +
                             ParameterBlockTypeName(Type()));
    try
    {
      return Value().GetValue<T>();
    }
    catch (const std::exception& exc)
    {
      throw std::logic_error(error_origin_scope_ + ":" + Name() + " " + exc.what());
    }
  }

  /**Fetches the parameter with the given name and returns it value.*/
  template <typename T>
  T GetParamValue(const std::string& param_name) const
  {
    try
    {
      const auto& param = GetParam(param_name);
      return param.GetValue<T>();
    }
    catch (const std::out_of_range& oor)
    {
      throw std::out_of_range(error_origin_scope_ + std::string(__PRETTY_FUNCTION__) +
                              ": Parameter \"" + param_name + "\" not present in block");
    }
  }

  /**Converts the parameters of an array-type parameter block to a vector of
   * primitive types and returns it.*/
  template <typename T>
  std::vector<T> GetVectorValue() const
  {
    if (Type() != ParameterBlockType::ARRAY)
      throw std::logic_error(error_origin_scope_ + std::string(__PRETTY_FUNCTION__) +
                             ": Invalid type requested for parameter of type " +
                             ParameterBlockTypeName(Type()));

    std::vector<T> vec;
    if (parameters_.empty()) return vec;

    // Check the first sub-param is of the right type
    const auto& front_param = parameters_.front();

    // Check that all other parameters are of the required type
    for (const auto& param : parameters_)
      if (param.Type() != front_param.Type())
        throw std::logic_error(error_origin_scope_ + " " + std::string(__PRETTY_FUNCTION__) +
                               ": Parameter \"" + name_ +
                               "\", cannot construct vector from block because "
                               "the sub_parameters do not all have the correct type. param->" +
                               ParameterBlockTypeName(param.Type()) + " vs param0->" +
                               ParameterBlockTypeName(front_param.Type()));

    const size_t num_params = parameters_.size();
    for (size_t k = 0; k < num_params; ++k)
    {
      const auto& param = GetParam(k);
      vec.push_back(param.GetValue<T>());
    }

    return vec;
  }

  /**Gets a vector of primitive types from an array-type parameter block
   * specified as a parameter of the current block.*/
  template <typename T>
  std::vector<T> GetParamVectorValue(const std::string& param_name) const
  {
    const auto& param = GetParam(param_name);
    return param.GetVectorValue<T>();
  }

  // Iterator
  // clang-format off
  class iterator
  {
  public:
    ParameterBlock& ref_block;
    size_t ref_id;

    iterator(ParameterBlock& in_block, size_t i) :
             ref_block(in_block), ref_id(i) {}

    iterator operator++()    { iterator i = *this; ref_id++; return i; }
    iterator operator++(int) { ref_id++; return *this; }

    ParameterBlock& operator*() { return ref_block.parameters_[ref_id]; }
    bool operator==(const iterator& rhs) const { return ref_id == rhs.ref_id; }
    bool operator!=(const iterator& rhs) const { return ref_id != rhs.ref_id; }
  };
  class const_iterator
  {
  public:
    const ParameterBlock& ref_block;
    size_t ref_id;

    const_iterator(const ParameterBlock& in_block, size_t i) :
                   ref_block(in_block), ref_id(i) {}

    const_iterator operator++(){ const_iterator i = *this; ref_id++; return i; }
    const_iterator operator++(int) { ref_id++; return *this; }

    const ParameterBlock& operator*() { return ref_block.parameters_[ref_id]; }
    bool operator==(const const_iterator& rhs) const
    { return ref_id == rhs.ref_id; }
    bool operator!=(const const_iterator& rhs) const
    { return ref_id != rhs.ref_id; }
  };
  // clang-format on

  iterator begin() { return {*this, 0}; }
  iterator end() { return {*this, parameters_.size()}; }

  const_iterator begin() const { return {*this, 0}; }
  const_iterator end() const { return {*this, parameters_.size()}; }

  /**Given a reference to a string, recursively travels the parameter
   * tree and print values into the reference string.*/
  void RecursiveDumpToString(std::string& outstr, const std::string& offset = "") const;

  void RecursiveDumpToJSON(std::string& outstr) const;
};

} // namespace chi
