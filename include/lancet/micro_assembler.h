#pragma once

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "absl/status/status.h"
#include "absl/time/time.h"
#include "absl/types/span.h"
#include "blockingconcurrentqueue.h"
#include "concurrentqueue.h"
#include "lancet/cli_params.h"
#include "lancet/read_extractor.h"
#include "lancet/ref_window.h"
#include "lancet/variant.h"
#include "lancet/variant_store.h"

namespace lancet {
struct WindowResult {
  absl::Duration runtime = absl::ZeroDuration();  // NOLINT
  std::size_t windowIdx = 0;                      // NOLINT

  [[nodiscard]] auto IsEmpty() const -> bool { return runtime == absl::ZeroDuration() && windowIdx == 0; }
};

using InWindowQueue = moodycamel::ConcurrentQueue<std::shared_ptr<RefWindow>>;
using OutResultQueue = moodycamel::BlockingConcurrentQueue<WindowResult>;

class MicroAssembler {
 public:
  explicit MicroAssembler(std::shared_ptr<InWindowQueue> winq, std::shared_ptr<OutResultQueue> resq,
                          std::shared_ptr<const CliParams> p)
      : windowQPtr(std::move(winq)), resultQPtr(std::move(resq)), params(std::move(p)) {}

  MicroAssembler() = default;

  void Process(const std::shared_ptr<VariantStore>& store);

 private:
  std::shared_ptr<InWindowQueue> windowQPtr;
  std::shared_ptr<OutResultQueue> resultQPtr;
  std::shared_ptr<const CliParams> params;

  std::vector<Variant> variants;
  std::vector<WindowResult> results;

  [[nodiscard]] auto ProcessWindow(ReadExtractor* re, const std::shared_ptr<const RefWindow>& w) -> absl::Status;
  [[nodiscard]] auto ShouldSkipWindow(const std::shared_ptr<const RefWindow>& w) const -> bool;

  // Try to flush variants to store if its possible to write to store without waiting for other threads
  void TryFlush(const std::shared_ptr<VariantStore>& store, const moodycamel::ProducerToken& token);

  // Force flush variants to store blocking current thread if any other threads are writing to store.
  void ForceFlush(const std::shared_ptr<VariantStore>& store, const moodycamel::ProducerToken& token);
};
}  // namespace lancet
