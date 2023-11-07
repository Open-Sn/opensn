#include "framework/math/chi_math.h"

#include "framework/math/Statistics/cdfsampler.h"

#include "framework/chi_runtime.h"
#include "framework/logging/chi_log.h"

#include <unistd.h>

/**Constructor for a sub interval*/
chi_math::CDFSampler::SubIntvl::SubIntvl(std::string offset,
                                         int ibin,
                                         int fbin,
                                         std::vector<double>& in_cdf,
                                         int subdiv_factor,
                                         int final_res,
                                         bool inhibit)
  : ref_cdf(in_cdf)
{
  inhibited = inhibit;
  cbin_i = ibin;
  cbin_f = fbin;

  if (!inhibited)
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

        if (i == (subdiv_factor - 1)) end = fbin;

        //        chi::log.Log()
        //          << offset
        //          << "Sub-interval " << beg
        //          << " " << end;

        sub_intvls[i] =
          new SubIntvl(offset + std::string("  "), beg, end, ref_cdf, subdiv_factor, final_res);
      }
    }
  }
}

/**Default constructor.*/
chi_math::CDFSampler::CDFSampler(std::vector<double>& in_cdf, int subdiv_factor, int final_res)
  : ref_cdf_(in_cdf)
{
  // Setting sub-division factor
  if (subdiv_factor >= 1) this->subdiv_factor_ = subdiv_factor;
  else
  {
    if (in_cdf.size() <= 10) this->subdiv_factor_ = 1;
    else if (in_cdf.size() <= 10000)
      this->subdiv_factor_ = 10;
    else
      this->subdiv_factor_ = 10; // sqrt(in_cdf.size());
  }

  // Setting final resolution
  if (final_res >= 3) this->final_res_ = final_res;
  else { this->final_res_ = 100; }

  //  chi::log.Log()
  //    << "Factors: " << this->subdiv_factor
  //    << " " << this->final_res;

  // Sub-dividing the interval
  size_t cdf_size = in_cdf.size();
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

      if (i == (this->subdiv_factor_ - 1)) end = cdf_size - 1;

      //      chi::log.Log()
      //        << "Sub-interval " << beg
      //        << " " << end;

      sub_intvls_[i] =
        new SubIntvl(std::string("  "), beg, end, ref_cdf_, this->subdiv_factor_, this->final_res_);
    }
  }
}

/**Initiates the sampling process.*/
int
chi_math::CDFSampler::Sample(double x)
{
  int ret_val = -1;
  int cdf_size = ref_cdf_.size();

  // Check bracket lo and hi
  if (x <= ref_cdf_[0]) ret_val = 0;
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
      if (sub_intvls_[s]->Sample(x, range)) break;
    }

    // Sample range
    //    chi::log.Log()
    //      << "Sampling " << x
    //      << " in range " << range.first << "-" << range.second;
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
    Chi::log.LogAllError() << "chi_math::CDFSampler::Sample. Error in CDF sampling routine. "
                           << "A bin was not found.";
    Chi::Exit(EXIT_FAILURE);
  }

  return ret_val;
}

/**Sampling a sub-interval.*/
bool
chi_math::CDFSampler::SubIntvl::Sample(double x, std::pair<int, int>& range)
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
      if (sub_intvls[s]->Sample(x, range)) return true;
    }
  }

  return false;
}

/**Sample a Cumulative Distribution Function (CDF) given a probability.
 *
 * The supplied vector should contain the upper bin boundary for each
 * bin and will return the bin associated with the bin that brackets
 * the supplied probability.
 *
 * Example:
 * Suppose we sample bins 0-9. Suppose also that the probalities for each
 * bin is as follows:
 * - 0.1 bin 0
 * - 0.1 bin 1
 * - 0.5 bin 5
 * - 0.3 bin 8
 *
 * The CDF for this probability distribution will look like this
 * - bin 0 = 0.1
 * - bin 1 = 0.2
 * - bin 2 = 0.2
 * - bin 3 = 0.2
 * - bin 4 = 0.2
 * - bin 5 = 0.7
 * - bin 6 = 0.7
 * - bin 7 = 0.7
 * - bin 8 = 1.0
 * - bin 9 = 1.0
 *
 * Supplying a random number between 0 and 1 should indicate sampling one
 * of the bins 0,1,5 or 8. The most inefficient way to do this is to
 * linearly loop through the cdf and check \f$ cdf_{i-1} \ge \theta < cdf_i \f$.
 *  An optimized version of this sampling would be to perform a recursive
 *  block search which starts with a course view of the cdf and then gradually
 *  refines the view until the final linear search can be performed.*/
int
chi_math::SampleCDF(double x, std::vector<double> cdf_bin)
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

  if ((indB - indA) <= fine_limit) refine_limit_reached = true;

  //  chi::log.Log() << "************ new prob x=" << x
  //   << " " << indA << " " << indB << " " << indC;
  //  refine_limit_reached = true;
  // Recursively refine
  int refine_count = 0;
  while (!refine_limit_reached)
  {
    int intvl_size = 0;
    if (x <= cdf_bin[indA]) refine_limit_reached = true;
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

    //    chi::log.Log()
    //      << refine_count << " "
    //      << "newA=" << indA
    //      << "newB=" << indB
    //      << "newC=" << indC;
    refine_count++;

    if (intvl_size <= fine_limit)
    {
      refine_limit_reached = true;
      lookup_i = indA;
      lookup_f = indC;
    }
  }

  //  chi::log.Log()
  //    << "Final lookup range=" << lookup_i << "-" << lookup_f << " " << x;

  // Perform final lookup
  int ret_val = -1;

  if (x <= cdf_bin[0]) ret_val = 0;
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
    Chi::log.LogAllError() << "chi_math::SampleCDF. Error in CDF sampling routine. "
                           << "A bin was not found."
                           << " i=" << lookup_i << " f=" << lookup_f << " x=" << x;
    Chi::Exit(EXIT_FAILURE);
  }

  //  chi::log.Log() << ret_val;
  //  usleep(100000);

  return ret_val;
}
