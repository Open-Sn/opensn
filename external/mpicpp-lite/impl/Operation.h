// SPDX-FileCopyrightText: 2023 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "mpi.h"

namespace mpicpp_lite {

namespace op {

/// Template for MPI operation `Op` on a `T` type
///
/// @tparam Op Operation
/// @tparam T Datatype
template <typename Op, typename T>
struct Operation {
};

// Sum

/// Template for summation operation on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct sum {
    /// Call operator
    ///
    /// @param x First operand
    /// @param y Second operand
    /// @return Sum of `x` and `y`
    const T &
    operator()(const T & x, const T & y) const
    {
        return x + y;
    }
};

/// Template for summation operation on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct Operation<sum<T>, T> {
    /// Call operator
    ///
    /// @return MPI operation for sumation
    static MPI_Op
    op()
    {
        return MPI_SUM;
    }
};

// Product

/// Template for product operation on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct prod {
    /// Call operator
    ///
    /// @param x First operand
    /// @param y Second operand
    /// @return Product of `x` and `y`
    const T &
    operator()(const T & x, const T & y) const
    {
        return x * y;
    }
};

/// Template for product operation on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct Operation<prod<T>, T> {
    /// Call operator
    ///
    /// @return MPI operation for product
    static MPI_Op
    op()
    {
        return MPI_PROD;
    }
};

// Maximum

/// Template for finding maximum on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct max {
    /// Call operator
    ///
    /// @param x First operand
    /// @param y Second operand
    /// @return The maximum of `x` and `y`
    const T &
    operator()(const T & x, const T & y) const
    {
        return x < y ? y : x;
    }
};

/// Template for finding maximum on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct Operation<max<T>, T> {
    /// Call operator
    ///
    /// @return MPI operation for finding maximum
    static MPI_Op
    op()
    {
        return MPI_MAX;
    }
};

// Minimum

/// Template for finding minimum on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct min {
    /// Call operator
    ///
    /// @param x First operand
    /// @param y Second operand
    /// @return The minimum of `x` and `y`
    const T &
    operator()(const T & x, const T & y) const
    {
        return x < y ? x : y;
    }
};

/// Template for finding minimum on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct Operation<min<T>, T> {
    /// Call operator
    ///
    /// @return MPI operation for finding minimum
    static MPI_Op
    op()
    {
        return MPI_MIN;
    }
};

// Logical AND

/// Template for logical AND on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct logical_and {
    /// Call operator
    ///
    /// @param x First operand
    /// @param y Second operand
    /// @return `x` AND `y`
    const T &
    operator()(const T & x, const T & y) const
    {
        return x && y;
    }
};

/// Template for logical AND on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct Operation<logical_and<T>, T> {
    /// Call operator
    ///
    /// @return MPI operation logical AND
    static MPI_Op
    op()
    {
        return MPI_LAND;
    }
};

// Logical OR

/// Template for logical OR on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct logical_or {
    /// Call operator
    ///
    /// @param x First operand
    /// @param y Second operand
    /// @return `x` OR `y`
    const T &
    operator()(const T & x, const T & y) const
    {
        return x || y;
    }
};

/// Template for logical OR on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct Operation<logical_or<T>, T> {
    /// Call operator
    ///
    /// @return MPI operation for logical OR
    static MPI_Op
    op()
    {
        return MPI_LOR;
    }
};

// Logical XOR

/// Template for logical XOR on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct logical_xor {
    /// Call operator
    ///
    /// @param x First operand
    /// @param y Second operand
    /// @return `x` XOR `y`
    const T &
    operator()(const T & x, const T & y) const
    {
        return !x != !y;
    }
};

/// Template for logical XOR on a `T` type
///
/// @tparam T Datatype
template <typename T>
struct Operation<logical_xor<T>, T> {
    /// Call operator
    ///
    /// @return MPI operation for logical XOR
    static MPI_Op
    op()
    {
        return MPI_LXOR;
    }
};

} // namespace op

} // namespace mpicpp_lite
