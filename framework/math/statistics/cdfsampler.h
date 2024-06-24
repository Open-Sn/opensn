// SPDX-FileCopyrightText: 2024 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

namespace opensn
{

/**Object for implementing an efficient cdf sampler.
 *
 * Normal linear sampling of bins of a Cumulative Distribution Function
 * is O(N) for each sample which can lead to very expensive sampling
 * of distributions with many samples. To this end this sampler
 * is designed to subdivide the bins recursively until a suitable
 * linear search can be performed.
 *
 * The speed-up for cdfs of >1000 bins is more than a factor 100.
 *
 * In order to use this sampler for repeated sampling calls make
 * sure to initialize it outside the phases that will repeatedly sample
 * it because it has some over-head to it that gets executed in the constructor.
 *
 \code
 CDFSampler sampler(cdf);
 \endcode
 *
 * */
class CDFSampler
{
public:
  struct SubIntvl;
  static const int AUTO_SUBDIV = -1;
  static const int AUTO_FINERES = -2;

private:
  int subdiv_factor_;
  int final_res_;
  std::vector<double>& ref_cdf_;
  std::vector<SubIntvl*> sub_intvls_;

public:
  /** constructor.*/
  CDFSampler(std::vector<double>& cdf,
             int subdiv_factor = AUTO_SUBDIV,
             int final_res = AUTO_FINERES);

  /**Initiates the sampling process.*/
  int Sample(double x);
};

// ###################################################################
/**Sub-structure for sub-intervals*/
struct CDFSampler::SubIntvl
{
  int cbin_i;
  int cbin_f;
  std::vector<double>& ref_cdf;
  bool inhibited;

  std::vector<std::shared_ptr<SubIntvl>> sub_intvls;

  std::string offset;

  /**Constructor for a sub interval*/
  SubIntvl(std::string offset,
           int ibin,
           int fbin,
           std::vector<double>& cdf,
           int subdiv_factor = 10,
           int final_res = 10,
           bool inhibit = false);

  /**Sampling a sub-interval.*/
  bool Sample(double x, std::pair<int, int>& range);
};

} // namespace opensn
