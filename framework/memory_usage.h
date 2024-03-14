#pragma once

namespace opensn
{

/**
 * Simple structure for memory usage.
 */
struct CSTMemory
{
  double memory_bytes = 0.0;
  double memory_kbytes = 0.0;
  double memory_mbytes = 0.0;
  double memory_gbytes = 0.0;

  CSTMemory() = default;

  explicit CSTMemory(double mem)
  {
    memory_bytes = mem;
    memory_kbytes = mem / 1024.0;
    memory_mbytes = mem / 1024.0 / 1024.0;
    memory_gbytes = mem / 1024.0 / 1024.0 / 1024.0;
  }

  CSTMemory& operator=(const CSTMemory& cst_mem) = default;
};

/**
 * Get current memory usage.
 */
CSTMemory GetMemoryUsage();

/**
 * Get current memory usage in Megabytes.
 */
double GetMemoryUsageInMB();

} // namespace opensn
