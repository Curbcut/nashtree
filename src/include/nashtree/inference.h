#ifndef NASHTREE_INFERENCE_H
#define NASHTREE_INFERENCE_H
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace nashtree {
// Explicit local-filesystem durability boundary (POSIX files and directories).
void sync_path(const std::string& path);
// Native core owns no R objects and performs no RNG or preprocessing operations.
struct MatrixView {
  const double* data;
  std::size_t rows, columns, row_stride, column_stride, buffer_size;
};
struct Node {
  std::int32_t feature, left, right;
  std::uint32_t kind; // 0 scalar leaf, 1 numeric split; NaN routes left.
  double threshold, value;
};
struct LoadLimits {
  std::uint64_t file_bytes = 2147483648ULL;
  std::uint64_t metadata_bytes = 16777216ULL;
  std::uint64_t nodes = 64000000ULL;
};
class Forest {
 public:
  static Forest load(const std::string& path, LoadLimits limits = {});
  void save(const std::string& path) const; // refuses replacement; atomic commit
  void validate() const;
  // Fixed output is column major: rows x requested IDs, including repeats.
  void predict_fixed(MatrixView x, const std::vector<std::uint64_t>& ids,
                     double* output, std::size_t output_size) const;
  // Mapped operation evaluates one explicitly designated sample per input row.
  void predict_mapped(MatrixView x, const std::vector<std::uint64_t>& ids,
                      double* output, std::size_t output_size) const;
  std::uint64_t features = 0, original_samples = 0;
  std::vector<std::uint64_t> sample_ids, sample_trees, tree_nodes;
  std::vector<Node> nodes;
  std::string metadata;
 private:
  std::size_t sample_index(std::uint64_t id) const;
  double evaluate(MatrixView x, std::size_t row, std::size_t sample) const;
  void validate_view(MatrixView x) const;
};
} // namespace nashtree
#endif
