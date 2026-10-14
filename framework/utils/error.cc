// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/utils/error.h"
#include "mpicpp-lite/mpicpp-lite.h"
#include <petscsys.h>
#include <sstream>
#include <type_traits>

namespace opensn
{

namespace
{

template <typename Return, typename Arg0, typename Arg1, typename Arg2>
Arg2 ThirdArgument(Return (*)(Arg0, Arg1, Arg2));

using PetscErrorMessageString = std::remove_pointer_t<decltype(ThirdArgument(&PetscErrorMessage))>;

} // namespace

[[noreturn]] void
ThrowMPIError(int ierr, const char* expr, const char* file, int line)
{
  std::stringstream message;
  message << "MPI call failed at " << file << ":" << line << " in expression \"" << expr
          << "\" with error code " << ierr << ". " << mpicpp_lite::error_message(ierr);

  throw std::runtime_error(message.str());
}

[[noreturn]] void
ThrowPETScError(int ierr, const char* expr, const char* file, int line)
{
  const char* ierr_desc = nullptr;
  PetscErrorMessageString ierr_desc_specific = nullptr;
  PetscErrorMessage(static_cast<PetscErrorCode>(ierr), &ierr_desc, &ierr_desc_specific);

  std::stringstream ss;
  ss << "PETSc call failed at " << file << ":" << line << " in expression \"" << expr
     << "\" with error code " << ierr;
  if (ierr_desc != nullptr)
    ss << ". " << ierr_desc;
  if (ierr_desc_specific != nullptr)
    ss << " (" << ierr_desc_specific << ")";

  throw std::runtime_error(ss.str());
}

} // namespace opensn
