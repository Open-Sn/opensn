// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#if defined(__NVCC__) || defined(__HIPCC__)
#define __inline_host_dev__ inline __host__ __device__
#else
#define __inline_host_dev__ inline
#endif
