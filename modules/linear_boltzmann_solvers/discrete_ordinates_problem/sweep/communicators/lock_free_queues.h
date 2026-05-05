// SPDX-FileCopyrightText: 2026 The OpenSn Authors <https://open-sn.github.io/opensn/>
// SPDX-License-Identifier: MIT

#pragma once

#include <atomic>
#include <cstddef>
#include <new>
#include <thread>
#include <vector>

namespace opensn
{

/**
 * Bounded lock-free multi-producer, single-consumer ring buffer.
 *
 * Producers reserve slots through an atomic head counter and publish them with a per-slot
 * ready flag. The single consumer drains in FIFO order through the tail index. The queue
 * is bounded and reuses preallocated slots; it performs no dynamic allocation once the
 * storage has been initialized.
 *
 * In the CBCD aggregated communicator, LockFreeRingBuffer serves two roles:
 * 1. an outgoing per-destination queue written by sweep worker threads and drained by the
 *    communication thread,
 * 2. an incoming per-angle-set queue written by the communication thread and drained by the
 *    owning angleset worker thread.
 *
 * LockFreeRingBuffer works under the following assumptions:
 * - producers reserve one slot, write the payload in place, and publish the slot exactly
 *   once
 * - the consumer drains published slots in FIFO order and returns them to the ring for
 *   reuse.
 *
 * This yields a fixed-capacity queue with explicit slot reuse.
 */
template <typename T>
class LockFreeRingBuffer
{
public:
  /// Slot payload with a publication flag.
  struct Slot
  {
    /// Stored payload.
    T payload;
    /// Publication flag visible to the single consumer.
    std::atomic<bool> ready{false};
  };

  /**
   * Allocate storage for the requested number of slots.
   *
   * \param capacity Number of ring-buffer slots.
   */
  void Preallocate(const std::size_t capacity) { buffer_ = std::vector<Slot>(capacity); }

  /**
   * Initialize every slot payload in place.
   *
   * \tparam Callback Callable invoked once per slot payload.
   * \param cb Initialization callback.
   */
  template <typename Callback>
  void InitializeSlots(Callback&& cb)
  {
    for (auto& slot : buffer_)
      cb(slot.payload);
  }

  /**
   * Reserve one slot for a producer.
   *
   * \return Writable slot reference.
   */
  Slot& ReserveSlot()
  {
    const auto idx = head_.fetch_add(1, std::memory_order_relaxed) % buffer_.size();
    while (buffer_[idx].ready.load(std::memory_order_acquire))
      std::this_thread::yield();
    return buffer_[idx];
  }

  /**
   * Publish one reserved slot to the consumer.
   *
   * \param slot Slot to publish.
   */
  void PublishSlot(Slot& slot) { slot.ready.store(true, std::memory_order_release); }

  /**
   * Gather currently ready slots without consuming them.
   *
   * \param out Output vector of ready slot pointers.
   */
  void GetReadySlots(std::vector<Slot*>& out)
  {
    out.clear();
    if (buffer_.empty())
      return;

    const auto capacity = buffer_.size();
    auto current_tail = tail_;
    while (buffer_[current_tail % capacity].ready.load(std::memory_order_acquire))
    {
      out.push_back(&buffer_[current_tail % capacity]);
      ++current_tail;
    }
  }

  /**
   * Release the next `count` ready slots after they have been consumed.
   *
   * \param count Number of slots to free.
   */
  void FreeSlots(const std::size_t count)
  {
    const auto capacity = buffer_.size();
    for (std::size_t i = 0; i < count; ++i)
    {
      buffer_[tail_ % capacity].ready.store(false, std::memory_order_release);
      ++tail_;
    }
  }

  /**
   * Consume all ready slots in FIFO order.
   *
   * \tparam Callback Callable invoked with each slot payload.
   * \param cb Consumer callback.
   * \return Number of consumed slots.
   */
  template <typename Callback>
  std::size_t ProcessReady(Callback&& cb)
  {
    if (buffer_.empty())
      return 0;

    const auto capacity = buffer_.size();
    std::size_t count = 0;
    while (true)
    {
      auto& slot = buffer_[tail_ % capacity];
      if (not slot.ready.load(std::memory_order_acquire))
        break;
      cb(slot.payload);
      slot.ready.store(false, std::memory_order_release);
      ++tail_;
      ++count;
    }
    return count;
  }

  /// Check whether the queue currently has no published slots.
  bool Empty() const
  {
    if (buffer_.empty())
      return true;
    return not buffer_[tail_ % buffer_.size()].ready.load(std::memory_order_acquire);
  }

private:
  /// Ring-buffer storage.
  std::vector<Slot> buffer_;
  /// Producer reservation index.
  alignas(std::hardware_destructive_interference_size) std::atomic<std::size_t> head_{0};
  /// Consumer drain index.
  alignas(std::hardware_destructive_interference_size) std::size_t tail_{0};
};

} // namespace opensn
