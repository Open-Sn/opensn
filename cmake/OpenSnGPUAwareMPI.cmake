# SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
# SPDX-License-Identifier: MIT

# Check that the MPI library can pass GPU buffers when OPENSN_USE_GPU_AWARE_MPI is ON.
#
# There is no portable way to ask an MPI library whether it is GPU-aware, so the check uses what
# each common implementation exposes:
#
# * Cray MPICH needs its GPU Transport Layer (GTL) library linked into the executable. Without it,
#   MPICH_GPU_SUPPORT_ENABLED=1 aborts in MPI_Init. The library is found and linked here.
# * MPICH reports its configure options through MPI_Get_library_version. A library configured
#   with --without-cuda, --without-hip, or --without-ze cannot pass GPU buffers. One configured
#   with --with-cuda, --with-hip, or --with-ze can. With neither option, MPICH enables a backend
#   only if its configure step detects the GPU runtime, and silently builds without it otherwise,
#   so support is reported as unverified.
# * Open MPI defines MPIX_CUDA_AWARE_SUPPORT and MPIX_ROCM_AWARE_SUPPORT in mpi-ext.h as 1 when it
#   was built with that support and as 0 when it was not. A library defining the macro as 0 is
#   rejected. One that does not define it is accepted with a warning.
#
# Other libraries are accepted with a message that GPU support could not be verified.
#
# The caller sets MPI_DETERMINE_LIBRARY_VERSION before find_package(MPI), so that FindMPI provides
# MPI_CXX_LIBRARY_VERSION_STRING. The check sets OPENSN_GPU_AWARE_MPI_LINK_LIBRARIES to the
# libraries the OpenSn library must link.

include(CheckCXXSourceCompiles)

function(opensn_check_gpu_aware_mpi)
  set(OPENSN_GPU_AWARE_MPI_LINK_LIBRARIES "" PARENT_SCOPE)

  if(OPENSN_WITH_CUDA OR HIP_ON_NVIDIA)
    set(backend "CUDA")
    set(gtl_name "mpi_gtl_cuda")
    set(mpich_backend "cuda")
    set(ompi_macro "MPIX_CUDA_AWARE_SUPPORT")
  elseif(OPENSN_WITH_HIP)
    set(backend "HIP")
    set(gtl_name "mpi_gtl_hsa")
    set(mpich_backend "hip")
    set(ompi_macro "MPIX_ROCM_AWARE_SUPPORT")
  else()
    set(backend "SYCL")
    set(gtl_name "mpi_gtl_ze")
    set(mpich_backend "ze")
    set(ompi_macro "")
  endif()

  set(version "${MPI_CXX_LIBRARY_VERSION_STRING}")
  if(NOT version OR version STREQUAL "NOTFOUND")
    set(version "${MPI_C_LIBRARY_VERSION_STRING}")
  endif()
  if(version STREQUAL "NOTFOUND")
    set(version "")
  endif()

  # Cray MPICH: link the GTL library for the GPU backend. The Cray environment variables identify
  # it only when the library version is unknown, since they stay set when another MPI is selected.
  if(version MATCHES "CRAY MPICH" OR
     (version STREQUAL "" AND (DEFINED ENV{CRAY_MPICH_ROOTDIR} OR DEFINED ENV{CRAY_MPICH_DIR})))
    # Cray installs the GTL in gtl/lib under the MPICH version directory, a few levels above the
    # MPI library (for example .../mpich/8.1.31/ofi/amd/6.0/lib/libmpi_amd.so). Site packages may
    # put it next to the MPI library instead, or link the MPI library from the Cray install. So
    # search the MPI library directory and gtl/lib above it, for the path FindMPI found and for
    # its resolved path.
    set(hints "")
    foreach(lib IN LISTS MPI_CXX_LIBRARIES MPI_C_LIBRARIES)
      if(NOT IS_ABSOLUTE "${lib}")
        continue()
      endif()
      get_filename_component(real "${lib}" REALPATH)
      foreach(path IN ITEMS "${lib}" "${real}")
        get_filename_component(dir "${path}" DIRECTORY)
        list(APPEND hints "${dir}")
        foreach(level RANGE 1 5)
          get_filename_component(dir "${dir}" DIRECTORY)
          list(APPEND hints "${dir}/gtl/lib")
        endforeach()
      endforeach()
    endforeach()
    foreach(var CRAY_MPICH_ROOTDIR CRAY_MPICH_DIR CRAY_MPICH_PREFIX)
      if(DEFINED ENV{${var}})
        list(APPEND hints "$ENV{${var}}/gtl/lib" "$ENV{${var}}/../../../gtl/lib")
      endif()
    endforeach()
    list(REMOVE_DUPLICATES hints)
    find_library(OPENSN_MPI_GTL_LIBRARY
      NAMES ${gtl_name}
      HINTS ${hints}
      DOC "Cray MPICH GPU Transport Layer library for GPU-aware MPI")
    if(NOT OPENSN_MPI_GTL_LIBRARY)
      list(JOIN hints "\n  " searched)
      message(FATAL_ERROR
        "OPENSN_USE_GPU_AWARE_MPI is ON with Cray MPICH, but lib${gtl_name} was not found. It is "
        "normally in gtl/lib under the Cray MPICH version directory. Set OPENSN_MPI_GTL_LIBRARY "
        "to its full path, or turn OPENSN_USE_GPU_AWARE_MPI OFF. Searched:\n  ${searched}")
    endif()
    # Nothing in OpenSn calls the GTL directly, so keep the linker from dropping it.
    set(OPENSN_GPU_AWARE_MPI_LINK_LIBRARIES
      "-Wl,--no-as-needed" "${OPENSN_MPI_GTL_LIBRARY}" "-Wl,--as-needed" PARENT_SCOPE)
    message(STATUS "GPU-aware MPI: Cray MPICH, linking ${OPENSN_MPI_GTL_LIBRARY}. Set "
                   "MPICH_GPU_SUPPORT_ENABLED=1 when running.")
    return()
  endif()

  # MPICH: the configure options show whether the backend was requested or disabled.
  if(version MATCHES "MPICH Version")
    set(end "([ \t\r\n'\"]|$)")
    if(version MATCHES "--without-${mpich_backend}${end}" OR
       version MATCHES "--with-${mpich_backend}=no${end}")
      message(FATAL_ERROR
        "OPENSN_USE_GPU_AWARE_MPI is ON, but the MPICH library was configured without "
        "${backend} support, so it cannot pass ${backend} buffers. Use an MPICH built with "
        "${backend} support, or turn OPENSN_USE_GPU_AWARE_MPI OFF.")
    endif()
    if(version MATCHES "--with-${mpich_backend}(=[^ \t\r\n'\"]*)?${end}")
      message(STATUS "GPU-aware MPI: MPICH was configured with --with-${mpich_backend}.")
    else()
      message(STATUS "GPU-aware MPI: could not verify that MPICH supports ${backend} buffers. Its "
                     "configure options do not include --with-${mpich_backend}, so ${backend} "
                     "support depends on what its configure step detected. Make sure it was built "
                     "with that support.")
    endif()
    return()
  endif()

  # Open MPI: mpi-ext.h defines the awareness macro as 1 when the support was built and as 0 when it
  # was not. Versions without the extension do not define it.
  if(version MATCHES "Open MPI" AND ompi_macro)
    set(CMAKE_REQUIRED_INCLUDES ${MPI_CXX_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${MPI_CXX_LIBRARIES})
    set(CMAKE_REQUIRED_QUIET ON)
    check_cxx_source_compiles("
      #include <mpi.h>
      #include <mpi-ext.h>
      #if !defined(${ompi_macro}) || !${ompi_macro}
      #error not ${backend}-aware
      #endif
      int main() { return 0; }" OPENSN_OPEN_MPI_${backend}_AWARE)
    if(OPENSN_OPEN_MPI_${backend}_AWARE)
      message(STATUS "GPU-aware MPI: Open MPI reports ${backend} support (${ompi_macro}).")
      return()
    endif()
    check_cxx_source_compiles("
      #include <mpi.h>
      #include <mpi-ext.h>
      #if !defined(${ompi_macro}) || ${ompi_macro}
      #error ${ompi_macro} is not defined as 0
      #endif
      int main() { return 0; }" OPENSN_OPEN_MPI_${backend}_UNAWARE)
    if(OPENSN_OPEN_MPI_${backend}_UNAWARE)
      message(FATAL_ERROR
        "OPENSN_USE_GPU_AWARE_MPI is ON, but Open MPI defines ${ompi_macro} as 0 in mpi-ext.h, so "
        "it was built without ${backend} support and cannot pass ${backend} buffers. Use an Open "
        "MPI built with ${backend} support, or turn OPENSN_USE_GPU_AWARE_MPI OFF.")
    endif()
    message(WARNING
      "OPENSN_USE_GPU_AWARE_MPI is ON, but ${ompi_macro} could not be read from Open MPI's "
      "mpi-ext.h, so ${backend} support could not be verified. If Open MPI was built without it, "
      "GPU-aware sweeps fail at run time.")
    return()
  endif()

  message(STATUS "GPU-aware MPI: could not verify that the MPI library supports ${backend} "
                 "buffers. Make sure it was built with that support.")
endfunction()
