#include "lancet/kmer.h"

#include "absl/hash/internal/city.h"
#include "lancet/merge_node_info.h"
#include "lancet/utils.h"

namespace lancet {
Kmer::Kmer(std::string_view sv) noexcept { Canonicalize(sv); }

auto Kmer::CanMergeKmers(const Kmer& buddy, BuddyPosition merge_dir, bool reverse_buddy, std::size_t k) const -> bool {
  return CanMergeSeqs(seq, buddy.SeqView(), merge_dir, reverse_buddy, k);
}

void Kmer::MergeBuddy(const Kmer& buddy, BuddyPosition dir, bool reverse_buddy, std::size_t k) {
  seq.reserve(seq.length() + buddy.Length() - k + 1);
  MergeKmerSeqs(&seq, buddy.seq, dir, reverse_buddy, k);
}

auto Kmer::FwdSeq() const -> std::string { return strand == Strand::REV ? utils::RevComp(seq) : seq; }

auto Kmer::IsCanonical(std::string_view sv) -> bool { return sv < utils::RevComp(sv); }

void Kmer::Canonicalize(std::string_view sv) {
  const auto revComp = utils::RevComp(sv);
  if (sv < revComp) {
    seq = std::string(sv);
    strand = Strand::FWD;
    return;
  }

  seq = revComp;
  strand = Strand::REV;
}

auto Kmer::ID() const -> std::uint64_t {
  return absl::hash_internal::CityHash64WithSeeds(seq.c_str(), seq.length(), utils::PRIME_0, utils::PRIME_1);  // NOLINT
}
}  // namespace lancet
