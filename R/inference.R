# Experimental raw-forest interface; deliberately not integrated into NASH runtime.
# Callers own preprocessing, feature order, transformations, random effects and RNG.
# IDs are original zero-based posterior sample IDs, not compact positions.
ntf_open <- function(path) inference_open_cpp(normalizePath(path, mustWork = TRUE))
ntf_close <- function(handle) invisible(inference_close_cpp(handle))
ntf_metadata <- function(handle) inference_metadata_cpp(handle)
ntf_matrix <- function(x) {
  stopifnot(is.matrix(x), is.numeric(x))
  storage.mode(x) <- "double"
  x
}
ntf_predict_fixed <- function(handle, x, sample_ids) {
  inference_fixed_cpp(handle, ntf_matrix(x), as.double(sample_ids))
}
ntf_predict_mapped <- function(handle, x, sample_ids) {
  inference_mapped_cpp(handle, ntf_matrix(x), as.double(sample_ids))
}

# Producer-neutral conversion boundary; all offsets and sample IDs are zero-based.
# Nodes are rows of feature,left,right,kind,threshold,value, in reachable preorder.
ntf_write <- function(path, features, original_samples, sample_ids,
                      sample_trees, tree_nodes, nodes, metadata) {
  stopifnot(is.character(path), length(path) == 1L, !is.na(path),
    is.character(metadata), length(metadata) == 1L, !is.na(metadata),
    is.matrix(nodes), is.numeric(nodes))
  inference_write_cpp(path, as.double(features), as.double(original_samples),
    as.double(sample_ids), as.double(sample_trees), as.double(tree_nodes),
    matrix(as.double(nodes), nrow(nodes), ncol(nodes)), metadata)
  invisible(path)
}

ntf_sync <- function(path) {
  stopifnot(is.character(path), length(path) == 1L, !is.na(path))
  inference_sync_cpp(path)
  invisible(path)
}
