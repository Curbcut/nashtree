# Build the inference-only R package from the SAME reviewed core/bindings.
# No training code, S3 registration, RNG, or duplicate evaluator is included.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(1L, 2L), !file.exists(args[[1]]))
revision <- if (length(args) == 2L) args[[2]] else system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
stopifnot(length(revision) == 1L, grepl("^[0-9a-f]{40}$", revision))
root <- normalizePath(".", mustWork = TRUE)
stopifnot(file.exists(file.path(root, "src/inference.cpp")))
out <- args[[1]]
dir.create(file.path(out, "src/include/nashtree"), recursive = TRUE)
dir.create(file.path(out, "R")); dir.create(file.path(out, "inst"))
files <- c("src/inference.cpp", "src/R_inference.cpp", "src/include/nashtree/inference.h", "R/inference.R")
for (p in files) stopifnot(file.copy(file.path(root, p), file.path(out, p)))
writeLines(c("Package: nashtreeInference", "Version: 0.1.2", "Title: Native Nashtree Inference Bindings",
  'Authors@R: person("David", "Wachsmuth", role = c("aut", "cre"), email = "david.wachsmuth@mcgill.ca")',
  "Description: Inference-only bindings to the nashtree_forest_v1 core. No training or S3 registrations.",
  "License: MIT + file LICENSE", "Encoding: UTF-8", "LinkingTo: cpp11", "Imports: digest, jsonlite", "SystemRequirements: C++17", paste0("Config/nashtree/SourceRevision: ", revision)), file.path(out, "DESCRIPTION"))
stopifnot(file.copy(file.path(root, "LICENSE"), file.path(out, "LICENSE")))
writeLines(c("CXX_STD = CXX17", "PKG_CXXFLAGS = -ffp-contract=off -fno-fast-math -ffile-prefix-map=$(CURDIR)=/nashtreeInference/src -fdebug-prefix-map=$(CURDIR)=/nashtreeInference/src"), file.path(out, "src/Makevars"))
writeLines(c('useDynLib(nashtreeInference, .registration = TRUE)',
  'export(ntf_open, ntf_close, ntf_metadata, ntf_predict_fixed, ntf_predict_mapped, ntf_write, ntf_sync)'), file.path(out, "NAMESPACE"))
writeLines(vapply(c(files, "src/cpp11.cpp", "R/cpp11.R", "native/r-package/build.R", "native/r-package/record-build.R"), function(p) paste(digest::digest(file = p, algo = "sha256"), p), character(1)),
  file.path(out, "inst/native-source-sha256.txt"))
# Subset the already generated registration, avoiding a second binding generator.
cpp <- readLines("src/cpp11.cpp")
starts <- which(cpp == "// R_inference.cpp")
stopifnot(length(starts) == 7L)
blocks <- unlist(lapply(starts, function(i) {
  end <- i + which(cpp[(i + 1L):length(cpp)] == "}")[[1L]]
  cpp[i:end]
}))
entries <- grep('_nashtree_inference_', cpp[grepl('DL_FUNC', cpp, fixed = TRUE)], value = TRUE, fixed = TRUE)
stopifnot(length(entries) == 7L)
writeLines(gsub("_nashtree_", "_nashtreeInference_", c(
  '#include "include/nashtree/inference.h"', '#include "cpp11/declarations.hpp"',
  '#include <R_ext/Visibility.h>', blocks,
  'static const R_CallMethodDef CallEntries[] = {', entries, '{NULL, NULL, 0}};',
  'extern "C" attribute_visible void R_init_nashtreeInference(DllInfo* dll) {',
  'R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);',
  'R_useDynamicSymbols(dll, FALSE); R_forceSymbols(dll, TRUE); }'), fixed = TRUE),
  file.path(out, "src/cpp11.cpp"))
r <- readLines("R/cpp11.R")
starts <- grep("^inference_.* <- function", r)
stopifnot(length(starts) == 7L)
writeLines(gsub("_nashtree_", "_nashtreeInference_", unlist(lapply(starts, function(i) {
  end <- i + which(r[(i + 1L):length(r)] == "}")[[1L]]
  r[i:end]
})), fixed = TRUE), file.path(out, "R/cpp11.R"))
cat(normalizePath(out), "\n")

# Capture actual compilation inputs before installation; retain with the package.
dir.create(file.path(out, "tools"))
stopifnot(file.copy("native/r-package/record-build.R", file.path(out, "tools/record-build.R")))
writeLines(c("#!/bin/sh", "set -eu", '\"${R_HOME}/bin/Rscript\" tools/record-build.R'), file.path(out, "configure"))
Sys.chmod(file.path(out, "configure"), mode = "0755")
