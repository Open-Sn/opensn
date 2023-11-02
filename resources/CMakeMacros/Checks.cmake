if(__CHECKS_INCLUDED)
    return()
endif()
set(__CHECKS_INCLUDED TRUE)

include(CheckTypeSize)
include(CheckSymbolExists)

CHECK_SYMBOL_EXISTS(PETSC_USE_64BIT_INDICES
                    "${PETSC_ROOT}/include/petscconf.h"
                    PETSC_USE_64BIT_INDICES)
if (NOT ${PETSC_USE_64BIT_INDICES} MATCHES 1)
    message(FATAL_ERROR "PETSc has not been configured with the flag "
                        " --with-64-bit-indices\n")
endif()

