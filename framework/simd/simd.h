// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>

#ifdef __AVX512F__
#include "framework/simd/simd_avx512.h"
#elif defined(__AVX__) || defined(__AVX2__)
#include "framework/simd/simd_avx.h"
#elif defined(__SSE2__) || defined(_M_X64) || defined(_M_AMD64) ||                                 \
  (defined(_M_IX86_FP) && (_M_IX86_FP >= 2))
#include "framework/simd/simd_sse2.h"
#elif (defined(__ARM_NEON) && (!defined(__ARM_ARCH) || (__ARM_ARCH >= 7))) ||                      \
  defined(__aarch64__) || defined(_M_ARM64) || defined(_M_ARM64EC)
#include "framework/simd/simd_neon.h"
#else
#include "framework/simd/simd_scalar.h"
#endif
