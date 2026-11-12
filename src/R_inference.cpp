#include <cpp11.hpp>
#include "include/nashtree/inference.h"
#include <cmath>
#include <climits>
#include <stdexcept>
namespace {
std::vector<std::uint64_t> native_ids(cpp11::doubles ids) {
  std::vector<std::uint64_t> out;out.reserve(ids.size());
  for(auto id:ids) {
    if(!std::isfinite(id)||id<0||id>9007199254740991.0||std::floor(id)!=id) cpp11::stop("invalid sample ID");
    out.push_back(static_cast<std::uint64_t>(id));
  }
  return out;
}
void require_open(cpp11::external_pointer<nashtree::Forest> ptr) {
  if(!ptr.get()) cpp11::stop("forest handle is closed");
}
}
[[cpp11::register]]
cpp11::external_pointer<nashtree::Forest> inference_open_cpp(std::string path) {
  return cpp11::external_pointer<nashtree::Forest>(new nashtree::Forest(nashtree::Forest::load(path)));
}
[[cpp11::register]]
void inference_close_cpp(cpp11::external_pointer<nashtree::Forest> ptr) { ptr.reset(); }
[[cpp11::register]]
std::string inference_metadata_cpp(cpp11::external_pointer<nashtree::Forest> ptr) {
  require_open(ptr);return ptr->metadata;
}
[[cpp11::register]]
cpp11::writable::doubles_matrix<> inference_fixed_cpp(cpp11::external_pointer<nashtree::Forest> ptr,
  cpp11::doubles_matrix<> x,cpp11::doubles ids) {
  require_open(ptr);auto selected=native_ids(ids);
  if(ids.size()>INT_MAX) cpp11::stop("too many selected samples");
  cpp11::writable::doubles_matrix<> out(x.nrow(),ids.size());
  ptr->predict_fixed({REAL(x),std::size_t(x.nrow()),std::size_t(x.ncol()),1,std::size_t(x.nrow()),std::size_t(x.size())},
    selected,REAL(out),out.size());return out;
}
[[cpp11::register]]
cpp11::writable::doubles inference_mapped_cpp(cpp11::external_pointer<nashtree::Forest> ptr,
  cpp11::doubles_matrix<> x,cpp11::doubles ids) {
  require_open(ptr);auto selected=native_ids(ids);cpp11::writable::doubles out(x.nrow());
  ptr->predict_mapped({REAL(x),std::size_t(x.nrow()),std::size_t(x.ncol()),1,std::size_t(x.nrow()),std::size_t(x.size())},
    selected,REAL(out),out.size());return out;
}

// A producer-neutral writer boundary: plain validated arrays only. A loaded
// training forest is traversed by its own package, never by casting its pointer
// into this inference library's ABI.
[[cpp11::register]]
void inference_write_cpp(std::string path, double features, double original_samples,
  cpp11::doubles ids, cpp11::doubles sample_trees, cpp11::doubles tree_nodes,
  cpp11::doubles_matrix<> nodes, std::string metadata) {
  auto scalar=[](double x) {
    if(!std::isfinite(x)||x<1||x>9007199254740991.0||std::floor(x)!=x)
      cpp11::stop("invalid forest dimension");
    return static_cast<std::uint64_t>(x);
  };
  if(nodes.ncol()!=6 || nodes.nrow()>64000000 || metadata.size()>16777216)
    cpp11::stop("forest writer input bound");
  const auto bytes=128.0+8.0*(ids.size()+sample_trees.size()+tree_nodes.size())+
    32.0*nodes.nrow()+metadata.size();
  if(bytes>2147483648.0) cpp11::stop("forest writer byte bound");
  nashtree::Forest out;
  out.features=scalar(features);out.original_samples=scalar(original_samples);
  out.sample_ids=native_ids(ids);out.sample_trees=native_ids(sample_trees);
  out.tree_nodes=native_ids(tree_nodes);out.metadata=std::move(metadata);
  out.nodes.reserve(nodes.nrow());
  for(int i=0;i<nodes.nrow();++i) {
    std::int32_t fields[4];
    for(int j=0;j<4;++j) {
      double x=nodes(i,j);
      if(!std::isfinite(x)||x<INT32_MIN||x>INT32_MAX||std::floor(x)!=x)
        cpp11::stop("invalid node integer");
      fields[j]=static_cast<std::int32_t>(x);
    }
    if(fields[3]<0||fields[3]>1) cpp11::stop("unsupported node kind");
    out.nodes.push_back({fields[0],fields[1],fields[2],static_cast<std::uint32_t>(fields[3]),nodes(i,4),nodes(i,5)});
  }
  out.save(path);
}

[[cpp11::register]]
void inference_sync_cpp(std::string path) { nashtree::sync_path(path); }
