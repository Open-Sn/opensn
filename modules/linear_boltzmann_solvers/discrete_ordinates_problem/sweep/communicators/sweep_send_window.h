// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <stdexcept>

namespace opensn
{

/**
 * Per-rank limit on posted sweep sends that have not completed.
 *
 * A rendezvous send holds a slot in the network transmit queue until the receiver matches it. When
 * a rank posts more sends than that queue holds, MPI can block inside the send call, and two ranks
 * blocked that way toward each other deadlock. A sender asks the window for slots before it posts a
 * unit of sends and returns them when the sends complete, so the rank never has more sends posted
 * than the limit.
 *
 * The angle sets of one groupset share a window. Groupsets sweep one at a time and every sweep
 * completes its sends before it returns, so the window bounds all posted sweep sends of the rank.
 */
class SweepSendWindow
{
public:
  /**
   * Window that admits at most `max_outstanding_sends` posted sends.
   *
   * \param max_outstanding_sends Maximum number of posted sends that have not completed. Zero
   *        means no limit.
   * \throw std::invalid_argument If `max_outstanding_sends` is negative.
   */
  explicit SweepSendWindow(int max_outstanding_sends)
    : max_outstanding_sends_(max_outstanding_sends)
  {
    if (max_outstanding_sends_ < 0)
      throw std::invalid_argument("SweepSendWindow: The send limit must not be negative.");
  }

  /**
   * Reserve slots for one unit of sends.
   *
   * A unit is admitted when no sends are outstanding, even if it is larger than the limit, so every
   * unit is eventually posted.
   *
   * \param count Number of sends in the unit.
   * \return True if the slots were reserved and the unit may be posted.
   * \throw std::invalid_argument If `count` is not positive.
   */
  bool TryAcquire(int count)
  {
    if (count <= 0)
      throw std::invalid_argument("SweepSendWindow: A send unit must hold at least one send.");
    if (max_outstanding_sends_ > 0 and num_outstanding_ > 0 and
        count > max_outstanding_sends_ - num_outstanding_)
      return false;
    num_outstanding_ += count;
    return true;
  }

  /**
   * Return the slots of completed sends.
   *
   * \param count Number of sends that completed.
   * \throw std::logic_error If `count` is negative or exceeds the outstanding sends.
   */
  void Release(int count)
  {
    if (count < 0 or count > num_outstanding_)
      throw std::logic_error("SweepSendWindow: Released more sends than are outstanding.");
    num_outstanding_ -= count;
  }

  /// Return the number of sends that are posted and not completed.
  int GetNumOutstanding() const { return num_outstanding_; }

private:
  /// Maximum number of outstanding sends, or zero for no limit.
  const int max_outstanding_sends_;
  /// Number of reserved slots whose sends have not completed.
  int num_outstanding_ = 0;
};

} // namespace opensn
