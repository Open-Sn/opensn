// SPDX-FileCopyrightText: 2023 David Andrs <andrsd@gmail.com>
// SPDX-License-Identifier: MIT

#pragma once

#include "mpi.h"

namespace mpicpp_lite {

/// Create a new datatype for MPI communication
///
/// @tparam T Datatype
/// @return New `MPI_Datatype`
template <typename T>
inline MPI_Datatype
build_mpi_datatype()
{
    return MPI_DATATYPE_NULL;
}

/// General template to obtain an MPI_Datatype from a C++ type
///
/// @tparam T C++ data type
/// @return `MPI_Datatype` that is used in the MPI API
template <typename T>
inline MPI_Datatype
get_mpi_datatype()
{
    return MPI_DATATYPE_NULL;
}

template <>
inline MPI_Datatype
get_mpi_datatype<char>()
{
    return MPI_BYTE;
}

template <>
inline MPI_Datatype
get_mpi_datatype<short>()
{
    return MPI_SHORT;
}

template <>
inline MPI_Datatype
get_mpi_datatype<int>()
{
    return MPI_INT;
}

template <>
inline MPI_Datatype
get_mpi_datatype<long int>()
{
    return MPI_LONG;
}

template <>
inline MPI_Datatype
get_mpi_datatype<long long int>()
{
    return MPI_LONG_LONG;
}

template <>
inline MPI_Datatype
get_mpi_datatype<unsigned char>()
{
    return MPI_UNSIGNED_CHAR;
}

template <>
inline MPI_Datatype
get_mpi_datatype<unsigned short>()
{
    return MPI_UNSIGNED_SHORT;
}

template <>
inline MPI_Datatype
get_mpi_datatype<unsigned int>()
{
    return MPI_UNSIGNED;
}

template <>
inline MPI_Datatype
get_mpi_datatype<unsigned long int>()
{
    return MPI_UNSIGNED_LONG;
}

template <>
inline MPI_Datatype
get_mpi_datatype<unsigned long long int>()
{
    return MPI_UNSIGNED_LONG_LONG;
}

template <>
inline MPI_Datatype
get_mpi_datatype<float>()
{
    return MPI_FLOAT;
}

template <>
inline MPI_Datatype
get_mpi_datatype<double>()
{
    return MPI_DOUBLE;
}

template <>
inline MPI_Datatype
get_mpi_datatype<long double>()
{
    return MPI_LONG_DOUBLE;
}

template <>
inline MPI_Datatype
get_mpi_datatype<bool>()
{
    return MPI_CXX_BOOL;
}

#if __cplusplus >= 201703L

template <>
inline MPI_Datatype
get_mpi_datatype<std::byte>()
{
    return MPI_BYTE;
}

#endif

} // namespace mpicpp_lite
