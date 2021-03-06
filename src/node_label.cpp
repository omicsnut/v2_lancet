#include "lancet/node_label.h"

#include <algorithm>

#include "absl/types/span.h"
#include "lancet/merge_node_info.h"

namespace lancet {
NodeLabel::NodeLabel(std::size_t node_len) { bases.resize(node_len); }

void NodeLabel::MergeBuddy(const NodeLabel& buddy, BuddyPosition dir, bool reverse_buddy, std::size_t k) {
  MergeNodeInfo(&bases, absl::MakeConstSpan(buddy.bases), dir, reverse_buddy, k);
}

void NodeLabel::Push(KmerLabel label) {
  const auto setLabel = [&label](BaseLabel& base) { base.SetLabel(label, true); };
  std::for_each(bases.begin(), bases.end(), setLabel);
}

auto NodeLabel::LabelRatio(KmerLabel label) const -> double {
  const auto hasLabel = [&label](const BaseLabel& base) { return base.HasLabel(label); };
  const auto count = std::count_if(bases.cbegin(), bases.cend(), hasLabel);
  return static_cast<double>(count) / static_cast<double>(bases.size());
}

auto NodeLabel::HasLabel(KmerLabel label) const -> bool {
  const auto hasLabel = [&label](const BaseLabel& base) { return base.HasLabel(label); };
  return std::any_of(bases.cbegin(), bases.cend(), hasLabel);
}

auto NodeLabel::IsLabelOnly(KmerLabel label) const -> bool {
  const auto hasLabel = [&label](const BaseLabel& base) { return base.HasLabel(label); };
  return std::all_of(bases.cbegin(), bases.cend(), hasLabel);
}

auto NodeLabel::FillColor() const -> std::string {
  const auto hasRef = HasLabel(KmerLabel::REFERENCE);
  const auto hasTmr = HasLabel(KmerLabel::TUMOR);
  const auto hasNml = HasLabel(KmerLabel::NORMAL);

  if (hasRef && hasTmr && hasNml) return "lightblue";
  if (hasTmr && !hasNml) return "orangered";
  if (hasNml && !hasTmr) return hasRef ? "lightblue" : "royalblue";
  return "lightblue";
}
}  // namespace lancet
