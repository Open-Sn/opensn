#include "console/chi_console.h"
#include "chi_console_structs.h"

#if defined(__MACH__)
#include <mach/mach.h>
#endif
#include <iostream>
#include <fstream>
#include <unistd.h>

//###################################################################
/**Gets the current memory usage.*/
chi::CSTMemory
chi::Console::GetMemoryUsage()
{
  double mem = 0.0;
#if defined(__MACH__)
  struct mach_task_basic_info info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  long long int bytes;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count) != KERN_SUCCESS)
  {
    bytes = 0;
  }
  bytes = info.resident_size;
  mem = (double)bytes;
#else
  long long int llmem = 0;
  long long int rss = 0;

  std::string ignore;
  std::ifstream ifs("/proc/self/stat", std::ios_base::in);
  ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >>
    ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >>
    ignore >> ignore >> ignore >> ignore >> llmem >> rss;

  long long int page_size_bytes = sysconf(_SC_PAGE_SIZE);
  mem = rss * page_size_bytes;
  /*
  FILE* fp = NULL;
  if((fp = fopen( "/proc/self/statm", "r" )) == NULL)
    return 0;
  if(fscanf(fp, "%*s%*s%*s%*s%*s%lld", &llmem) != 1)
  {
    fclose(fp);
    return 0;
  }
  fclose(fp);*/

  // mem = llmem * (long long int)sysconf(_SC_PAGESIZE);
#endif

  CSTMemory mem_struct(mem);

  return mem_struct;
}

//###################################################################
/**Gets the current memory usage in megabytes.*/
double
chi::Console::GetMemoryUsageInMB()
{
  CSTMemory mem_struct = GetMemoryUsage();

  return mem_struct.memory_mbytes;
}
