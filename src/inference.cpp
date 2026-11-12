#include "include/nashtree/inference.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

namespace nashtree {
void sync_path(const std::string& path) {
#ifdef _WIN32
  (void)path;
  throw std::runtime_error("durable conversion cache requires a supported POSIX filesystem");
#else
  const int fd=::open(path.c_str(),O_RDONLY);
  if(fd<0) throw std::runtime_error("cannot open durability target");
  const int status=::fsync(fd);
  const int closed=::close(fd);
  if(status!=0 || closed!=0) throw std::runtime_error("cannot synchronize durability target");
#endif
}
namespace {
void require(bool ok, const char* message) {
  if (!ok) throw std::runtime_error(message);
}
std::uint64_t add(std::uint64_t a, std::uint64_t b) {
  require(b <= UINT64_MAX - a, "integer size overflow"); return a + b;
}
std::uint64_t mul(std::uint64_t a, std::uint64_t b) {
  require(a == 0 || b <= UINT64_MAX / a, "integer size overflow"); return a * b;
}
std::uint64_t crc64(const unsigned char* p, std::size_t n, std::uint64_t crc = 0) {
  static const auto table = [] {
    std::array<std::uint64_t, 256> t{};
    for (std::uint64_t i = 0; i < 256; ++i) {
      auto v = i << 56;
      for (int k = 0; k < 8; ++k) v = (v << 1) ^ ((v >> 63) ? 0x42F0E1EBA9EA3693ULL : 0);
      t[i] = v;
    }
    return t;
  }();
  for (std::size_t i = 0; i < n; ++i) crc = table[(crc >> 56) ^ p[i]] ^ (crc << 8);
  return crc;
}
void put(std::vector<unsigned char>& b, std::uint64_t v, unsigned width = 8) {
  for (unsigned i = 0; i < width; ++i) b.push_back(static_cast<unsigned char>(v >> (8 * i)));
}
std::uint64_t get(const unsigned char* p, unsigned width = 8) {
  std::uint64_t v = 0;
  for (unsigned i = 0; i < width; ++i) v |= std::uint64_t(p[i]) << (8 * i);
  return v;
}
void put_double(std::vector<unsigned char>& b, double x) {
  std::uint64_t v; std::memcpy(&v, &x, 8); put(b, v);
}
double get_double(const unsigned char* p) {
  auto v = get(p); double x; std::memcpy(&x, &v, 8); return x;
}
constexpr std::size_t header_bytes = 128;
const std::array<unsigned char, 8> magic{'N','A','S','H','N','T','F',0};
}

void Forest::validate() const {
  require(sizeof(double) == 8 && std::numeric_limits<double>::is_iec559, "binary64 platform required");
  require(features > 0 && features <= INT32_MAX && original_samples > 0, "invalid dimensions");
  require(!sample_ids.empty() && sample_trees.size() == sample_ids.size() + 1, "invalid sample offsets");
  require(!tree_nodes.empty() && sample_trees.front() == 0 && sample_trees.back() == tree_nodes.size() - 1,
          "invalid sample tree ranges");
  require(tree_nodes.front() == 0 && tree_nodes.back() == nodes.size(), "invalid tree node ranges");
  require(!metadata.empty(), "metadata required");
  for (std::size_t s = 0; s < sample_ids.size(); ++s) {
    require(sample_ids[s] < original_samples && (s == 0 || sample_ids[s-1] < sample_ids[s]), "sample IDs must be unique and sorted");
    require(sample_trees[s] < sample_trees[s+1], "empty or reversed sample");
  }
  for (std::size_t t = 0; t + 1 < tree_nodes.size(); ++t) {
    const auto begin = tree_nodes[t], end = tree_nodes[t+1];
    require(begin < end && end <= nodes.size() && end - begin <= INT32_MAX, "invalid tree offsets");
    // Canonical preorder requires each child after its parent. Single-parent and
    // reachability checks reject DAGs, cycles, and trailing unreachable nodes.
    std::vector<unsigned char> parents(end - begin, 0);
    for (auto i = begin; i < end; ++i) {
      const auto& n = nodes[i];
      if (n.kind == 0) {
        require(n.feature == -1 && n.left == -1 && n.right == -1 && std::isfinite(n.value), "invalid scalar leaf");
      } else {
        require(n.kind == 1 && n.feature >= 0 && std::uint64_t(n.feature) < features && std::isfinite(n.threshold), "unsupported split or threshold");
        for (auto child : {n.left, n.right}) {
          require(child > std::int64_t(i - begin) && std::uint64_t(child) < end - begin, "invalid child or noncanonical tree");
          require(++parents[child] == 1, "tree has multiple parents");
        }
      }
    }
    for (std::size_t i = 1; i < parents.size(); ++i) require(parents[i] == 1, "unreachable node");
  }
}

Forest Forest::load(const std::string& path, LoadLimits limits) {
  const auto size = std::filesystem::file_size(path);
  require(size >= header_bytes && size <= limits.file_bytes && size <= SIZE_MAX, "file outside load limits");
  std::ifstream file(path, std::ios::binary);
  std::array<unsigned char, header_bytes> h{};
  require(bool(file.read(reinterpret_cast<char*>(h.data()), h.size())), "truncated header");
  require(std::equal(magic.begin(), magic.end(), h.begin()), "invalid magic");
  require(get(h.data()+8,4) == 1 && get(h.data()+12,4) == 0 && get(h.data()+16,4) == 0x01020304 && get(h.data()+20,4) == 1,
          "unsupported version, byte order or capabilities");
  require(get(h.data()+24) == size, "file length mismatch");
  const auto ns = get(h.data()+32), nt = get(h.data()+40), nn = get(h.data()+48), nm = get(h.data()+56);
  require(ns > 0 && nt >= ns && nn >= nt && nn <= limits.nodes && nm > 0 && nm <= limits.metadata_bytes, "section counts outside load limits");
  const auto ids_offset = header_bytes;
  const auto samples_offset = add(ids_offset, mul(ns, 8));
  const auto trees_offset = add(samples_offset, mul(add(ns, 1), 8));
  const auto nodes_offset = add(trees_offset, mul(add(nt, 1), 8));
  const auto metadata_offset = add(nodes_offset, mul(nn, 32));
  require(add(metadata_offset, nm) == size && get(h.data()+88) == ids_offset && get(h.data()+96) == samples_offset &&
          get(h.data()+104) == trees_offset && get(h.data()+112) == nodes_offset && get(h.data()+120) == metadata_offset,
          "invalid section offsets");
  std::vector<unsigned char> body(size - header_bytes);
  require(bool(file.read(reinterpret_cast<char*>(body.data()), body.size())), "truncated payload");
  // CRC protects header, payload and provenance metadata (checksum slot zeroed).
  // The caller additionally verifies the immutable artifact SHA-256.
  const auto expected_crc=get(h.data()+80);
  std::fill(h.begin()+80,h.begin()+88,0);
  require(crc64(body.data(),body.size(),crc64(h.data(),h.size()))==expected_crc,"file checksum mismatch");
  Forest result; result.features = get(h.data()+64); result.original_samples = get(h.data()+72);
  auto ptr = [&](std::uint64_t at) { return body.data() + at - header_bytes; };
  result.sample_ids.reserve(ns); result.sample_trees.reserve(ns+1); result.tree_nodes.reserve(nt+1); result.nodes.reserve(nn);
  for (std::uint64_t i=0;i<ns;++i) result.sample_ids.push_back(get(ptr(ids_offset+8*i)));
  for (std::uint64_t i=0;i<=ns;++i) result.sample_trees.push_back(get(ptr(samples_offset+8*i)));
  for (std::uint64_t i=0;i<=nt;++i) result.tree_nodes.push_back(get(ptr(trees_offset+8*i)));
  for (std::uint64_t i=0;i<nn;++i) {
    auto p = ptr(nodes_offset+32*i);
    result.nodes.push_back({static_cast<std::int32_t>(get(p,4)), static_cast<std::int32_t>(get(p+4,4)),
      static_cast<std::int32_t>(get(p+8,4)), static_cast<std::uint32_t>(get(p+12,4)), get_double(p+16), get_double(p+24)});
  }
  result.metadata.assign(reinterpret_cast<const char*>(ptr(metadata_offset)), nm);
  result.validate(); return result;
}

void Forest::save(const std::string& path) const {
  validate();
  require(!std::filesystem::exists(path), "refusing to replace artifact");
  const auto ns=sample_ids.size(), nt=tree_nodes.size()-1, nn=nodes.size();
  std::vector<unsigned char> body;
  for (auto v : sample_ids) put(body,v);
  for (auto v : sample_trees) put(body,v);
  for (auto v : tree_nodes) put(body,v);
  for (const auto& n : nodes) {
    put(body,static_cast<std::uint32_t>(n.feature),4); put(body,static_cast<std::uint32_t>(n.left),4);
    put(body,static_cast<std::uint32_t>(n.right),4); put(body,n.kind,4); put_double(body,n.threshold); put_double(body,n.value);
  }
  body.insert(body.end(),metadata.begin(),metadata.end());
  std::vector<unsigned char> h(magic.begin(),magic.end());
  put(h,1,4); put(h,0,4); put(h,0x01020304,4); put(h,1,4);
  put(h,add(header_bytes,body.size())); put(h,ns); put(h,nt); put(h,nn); put(h,metadata.size());
  put(h,features); put(h,original_samples); put(h,0);
  put(h,header_bytes); put(h,header_bytes+8*ns); put(h,header_bytes+8*ns+8*(ns+1));
  put(h,header_bytes+8*ns+8*(ns+1)+8*(nt+1)); put(h,header_bytes+8*ns+8*(ns+1)+8*(nt+1)+32*nn);
  require(h.size()==header_bytes, "internal header size");
  const auto checksum=crc64(body.data(),body.size(),crc64(h.data(),h.size()));
  for(unsigned i=0;i<8;++i) h[80+i]=static_cast<unsigned char>(checksum>>(8*i));
  // Unique temporary sibling, then no-replace hard-link commit on same filesystem.
  const auto staging=path+".writing";
  require(std::filesystem::create_directory(staging), "output writer already active");
  const auto temp=(std::filesystem::path(staging)/"artifact").string();
  try {
    std::ofstream out(temp,std::ios::binary);
    require(bool(out.write(reinterpret_cast<const char*>(h.data()),h.size())) &&
            bool(out.write(reinterpret_cast<const char*>(body.data()),body.size())), "write failed");
    out.close(); require(bool(out), "close failed");
    // Validate bytes read from disk before exposing a completed artifact.
    Forest::load(temp);
    std::filesystem::create_hard_link(temp,path);
    std::filesystem::remove_all(staging);
  } catch (...) { std::filesystem::remove_all(staging); throw; }
}

std::size_t Forest::sample_index(std::uint64_t id) const {
  const auto it=std::lower_bound(sample_ids.begin(),sample_ids.end(),id);
  require(it!=sample_ids.end() && *it==id, "requested sample was not retained");
  return it-sample_ids.begin();
}
void Forest::validate_view(MatrixView x) const {
  require(x.columns==features, "feature dimension mismatch");
  if (x.rows==0) return;
  require(x.data && x.row_stride>0 && x.column_stride>0, "invalid input view");
  require(add(mul(x.rows-1,x.row_stride),mul(x.columns-1,x.column_stride))<x.buffer_size, "input view outside buffer");
}
double Forest::evaluate(MatrixView x, std::size_t row, std::size_t sample) const {
  double sum=0.0;
  for (auto t=sample_trees[sample];t<sample_trees[sample+1];++t) {
    const auto root=tree_nodes[t]; auto at=root;
    while (nodes[at].kind!=0) {
      const auto& n=nodes[at]; const auto v=x.data[row*x.row_stride+n.feature*x.column_stride];
      at=root+((std::isnan(v)||v<=n.threshold)?n.left:n.right);
    }
    sum+=nodes[at].value; // Exact legacy tree order; no reassociated reduction.
  }
  return sum;
}
void Forest::predict_fixed(MatrixView x,const std::vector<std::uint64_t>& ids,double* out,std::size_t size) const {
  validate_view(x); require(mul(x.rows,ids.size())==size && (size==0||out), "output buffer mismatch");
  std::vector<std::size_t> indices; for(auto id:ids) indices.push_back(sample_index(id));
  for(std::size_t s=0;s<ids.size();++s) for(std::size_t r=0;r<x.rows;++r) out[s*x.rows+r]=evaluate(x,r,indices[s]);
}
void Forest::predict_mapped(MatrixView x,const std::vector<std::uint64_t>& ids,double* out,std::size_t size) const {
  validate_view(x); require(ids.size()==x.rows && size==x.rows && (size==0||out), "mapped buffer mismatch");
  std::vector<std::size_t> indices; for(auto id:ids) indices.push_back(sample_index(id));
  for(std::size_t r=0;r<x.rows;++r) out[r]=evaluate(x,r,indices[r]);
}
} // namespace nashtree
