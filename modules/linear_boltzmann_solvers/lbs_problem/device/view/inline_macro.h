// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#if defined(__NVCC__) || defined(__HIPCC__)
#define OPENSN_INLINE_HOST_DEV inline __host__ __device__
#else
#define OPENSN_INLINE_HOST_DEV inline
#endif
