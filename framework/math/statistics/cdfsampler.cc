// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#include "framework/math/math.h"

#include "framework/math/statistics/cdfsampler.h"

#include "framework/runtime.h"
#include "framework/logging/log.h"

#include <unistd.h>

namespace opensn
{

CDFSampler::SubIntvl::SubIntvl(std::string offset,
                               int ibin,
                               int fbin,
                               std::vector<double>& cdf,
                               int subdiv_factor,
                               int final_res,
                               bool inhibit)
  : ref_cdf(cdf)
{
  inhibited = inhibit;
  cbin_i = ibin;
  cbin_f = fbin;

  if (not inhibited)
  {
    size_t cdf_size = cbin_f - cbin_i + 1;
    size_t intvl_size = ceil(cdf_size / (double)subdiv_factor);

    if (intvl_size < final_res)
    {
      sub_intvls.push_back(new SubIntvl(
        offset + std::string("  "), ibin, fbin, ref_cdf, subdiv_factor, final_res, true));
    }
    else
    {
      sub_intvls.resize(subdiv_factor);
      for (int i = 0; i < subdiv_factor; i++)
      {
        int beg = ibin + i * intvl_size;
        int end = ibin + (i + 1) * intvl_size - 1;

        if (i == (subdiv_factor - 1))
          end = fbin;

        sub_intvls[i] =
          new SubIntvl(offset + std::string("  "), beg, end, ref_cdf, subdiv_factor, final_res);
      }
    }
  }
}

CDFSampler::CDFSampler(std::vector<double>& cdf, int subdiv_factor, int final_res) : ref_cdf_(cdf)
{
  // Setting sub-division factor
  if (subdiv_factor >= 1)
    this->subdiv_factor_ = subdiv_factor;
  else
  {
    if (cdf.size() <= 10)
      this->subdiv_factor_ = 1;
    else if (cdf.size() <= 10000)
      this->subdiv_factor_ = 10;
    else
      this->subdiv_factor_ = 10; // sqrt(cdf.size());
  }

  // Setting final resolution
  if (final_res >= 3)
    this->final_res_ = final_res;
  else
  {
    this->final_res_ = 100;
  }

  // Sub-dividing the interval
  size_t cdf_size = cdf.size();
  size_t intvl_size = ceil(cdf_size / (double)this->subdiv_factor_);

  if (intvl_size < this->final_res_)
    sub_intvls_.push_back(new SubIntvl(
      std::string("  "), 0, cdf_size - 1, ref_cdf_, this->subdiv_factor_, this->final_res_, true));
  else
  {
    sub_intvls_.resize(this->subdiv_factor_);
    for (int i = 0; i < this->subdiv_factor_; i++)
    {
      int beg = i * intvl_size;
      int end = (i + 1) * intvl_size - 1;

      if (i == (this->subdiv_factor_ - 1))
        end = cdf_size - 1;

      sub_intvls_[i] =
        new SubIntvl(std::string("  "), beg, end, ref_cdf_, this->subdiv_factor_, this->final_res_);
    }
  }
}

int
CDFSampler::Sample(double x)
{
  int ret_val = -1;
  int cdf_size = ref_cdf_.size();

  // Check bracket lo and hi
  if (x <= ref_cdf_[0])
    ret_val = 0;
  else if (x >= ref_cdf_[cdf_size - 1])
    ret_val = cdf_size - 1;
  // Check internal
  else
  {
    std::pair<int, int> range(0, cdf_size - 1);

    // Sample sub-intvls for range
    int num_sub_intvls = sub_intvls_.size();
    for (int s = 0; s < num_sub_intvls; s++)
    {
      if (sub_intvls_[s]->Sample(x, range))
        break;
    }

    for (int k = range.first; k <= range.second; k++)
    {
      if (k == 0)
      {
        if (x < ref_cdf_[k])
        {
          ret_val = k;
          break;
        }
      }
      else if ((x >= ref_cdf_[k - 1]) and (x < ref_cdf_[k]))
      {
        ret_val = k;
        break;
      }
    } // for k

  } // if internal

  if (ret_val < 0)
  {
    log.LogAllError() << "CDFSampler::Sample. Error in CDF sampling routine. "
                      << "A bin was not found.";
    Exit(EXIT_FAILURE);
  }

  return ret_val;
}

bool
CDFSampler::SubIntvl::Sample(double x, std::pair<int, int>& range)
{
  // If this was an inhibited intvl
  if (inhibited)
  {
    if (cbin_i == 0)
    {
      if (x < ref_cdf[cbin_i])
      {
        range.first = cbin_i;
        range.second = cbin_f;

        return true;
      }
    }

    if ((x >= ref_cdf[cbin_i - 1]) and (x < ref_cdf[cbin_f]))
    {
      range.first = cbin_i;
      range.second = cbin_f;

      return true;
    }
  }
  // If not inhibited sample sub-intvls
  else
  {
    int num_sub_intvls = sub_intvls.size();
    for (int s = 0; s < num_sub_intvls; s++)
    {
      if (sub_intvls[s]->Sample(x, range))
        return true;
    }
  }

  return false;
}

int
SampleCDF(double x, std::vector<double> cdf_bin)
{
  size_t fine_limit = 5;
  size_t cdf_size = cdf_bin.size();

  size_t lookup_i = 0;
  size_t lookup_f = cdf_size - 1;

  // Initial coursest level
  size_t indA = 0;
  size_t indB = std::ceil(cdf_size / 2.0) - 1;
  size_t indC = cdf_size - 1;

  bool refine_limit_reached = false;

  if ((indB - indA) <= fine_limit)
    refine_limit_reached = true;

  // Recursively refine
  int refine_count = 0;
  while (not refine_limit_reached)
  {
    int intvl_size = 0;
    if (x <= cdf_bin[indA])
      refine_limit_reached = true;
    else if (x > cdf_bin[indC])
      refine_limit_reached = true;
    else if ((x >= cdf_bin[indA]) and (x < cdf_bin[indB]))
    {
      intvl_size = indB - indA + 1;

      indC = indB;
      indB = indA + std::ceil(intvl_size / 2.0) - 1;
    }
    else
    {
      intvl_size = indC - indB + 1;

      indA = indB;
      indB = indA + std::ceil(intvl_size / 2.0) - 1;
    }

    refine_count++;

    if (intvl_size <= fine_limit)
    {
      refine_limit_reached = true;
      lookup_i = indA;
      lookup_f = indC;
    }
  }

  // Perform final lookup
  int ret_val = -1;

  if (x <= cdf_bin[0])
    ret_val = 0;
  else if (x >= cdf_bin[cdf_size - 1])
    ret_val = cdf_size - 1;
  else
  {
    for (int k = lookup_i; k <= lookup_f; k++)
    {
      if (k == 0)
      {
        if (x < cdf_bin[k])
        {
          ret_val = k;
          break;
        }
      }
      else if ((x >= cdf_bin[k - 1]) and (x < cdf_bin[k]))
      {
        ret_val = k;
        break;
      }
    } // for k
  }

  if (ret_val < 0)
  {
    log.LogAllError() << "SampleCDF. Error in CDF sampling routine. "
                      << "A bin was not found."
                      << " i=" << lookup_i << " f=" << lookup_f << " x=" << x;
    Exit(EXIT_FAILURE);
  }

  return ret_val;
}

} // namespace opensn
