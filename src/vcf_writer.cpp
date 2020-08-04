#include "lancet/vcf_writer.h"

#include <stdexcept>

#include "absl/strings/str_format.h"
#include "htslib/bgzf.h"
#include "htslib/tbx.h"

namespace lancet {
class VcfWriter::Impl {
 public:
  explicit Impl(const std::filesystem::path& out_path) : vcfPath(out_path) {
    fp = bgzf_open(out_path.c_str(), "w");
    if (fp == nullptr) {
      const auto errMsg = absl::StrFormat("could not open BGZF handle for %s", out_path);
      throw std::runtime_error(errMsg);
    }
  }

  auto Write(std::string_view record) -> absl::Status {
    const auto numWritten = bgzf_write(fp, record.data(), record.length());
    return numWritten == record.length() ? absl::OkStatus() : absl::InternalError("could not write to BGZF handle");
  }

  void Close() {
    if (fp == nullptr || isClosed) return;
    bgzf_close(fp);
    tbx_index_build(vcfPath.c_str(), 0, &tbx_conf_vcf);
  }

 private:
  bool isClosed = false;
  std::filesystem::path vcfPath;
  BGZF* fp = nullptr;
};

VcfWriter::VcfWriter(const std::filesystem::path& out_path) : pimpl(std::make_unique<Impl>(out_path)) {}
VcfWriter::~VcfWriter() { Close(); }
auto VcfWriter::Write(std::string_view record) -> absl::Status { return pimpl->Write(record); }
void VcfWriter::Close() { return pimpl->Close(); }
}  // namespace lancet