#pragma once

//=============================================================================
namespace chi
{
/**Simple structure for memory usage.*/
struct CSTMemory
{
  double memory_bytes = 0.0;
  double memory_kbytes = 0.0;
  double memory_mbytes = 0.0;
  double memory_gbytes = 0.0;

  CSTMemory() = default;

  explicit CSTMemory(double in_mem)
  {
    memory_bytes = in_mem;
    memory_kbytes = in_mem / 1024.0;
    memory_mbytes = in_mem / 1024.0 / 1024.0;
    memory_gbytes = in_mem / 1024.0 / 1024.0;
  }

  CSTMemory& operator=(const CSTMemory& in_struct) = default;
};
} // namespace chi
