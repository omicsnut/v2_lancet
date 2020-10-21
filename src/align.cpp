#include "lancet/align.h"

#include <algorithm>
#include <utility>

#include "needleman_wunsch.h"

namespace lancet {

auto Align(std::string_view ref, std::string_view qry) -> AlignedSequences {
  static scoring_t scoring;
  scoring_system_default(&scoring);

  nw_aligner_t* aligner = needleman_wunsch_new();
  alignment_t* aln = alignment_create(std::max(ref.size(), qry.size()));
  needleman_wunsch_align2(ref.data(), qry.data(), ref.size(), qry.size(), &scoring, aligner, aln);

  std::string alnRef(aln->result_a);
  std::string alnQry(aln->result_b);

  alignment_free(aln);
  needleman_wunsch_free(aligner);
  return AlignedSequences(std::move(alnRef), std::move(alnQry));
}

auto TrimEndGaps(AlignedSequencesView* aln) -> std::size_t {
  // Trim end GAPS and adjust end alignments until both ends in ref and qry have no GAPS
  std::size_t refStartTrim = 0;
  std::size_t start = 0;
  std::size_t end = aln->ref.length() - 1;

  const auto startGap = aln->ref[start] == ALIGN_GAP || aln->qry[start] == ALIGN_GAP;
  const auto endGap = aln->ref[end] == ALIGN_GAP || aln->qry[end] == ALIGN_GAP;

  if (startGap || endGap) {
    // move start until no begin alignment gaps are found
    while (aln->ref[start] == ALIGN_GAP || aln->qry[start] == ALIGN_GAP) {
      if (aln->ref[start] == ALIGN_GAP) refStartTrim++;
      start++;
    }

    // move end until no end alignment gaps are found
    while (aln->ref[end] == ALIGN_GAP || aln->qry[end] == ALIGN_GAP) {
      end--;
    }

    aln->ref = aln->ref.substr(start, end - start);
    aln->qry = aln->qry.substr(start, end - start);
  }

  return refStartTrim;
}
}  // namespace lancet
