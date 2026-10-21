# configure runs in the package source directory, before compilation.
config <- function(key) {
  x <- system2(file.path(R.home("bin"), "R"), c("CMD", "config", key), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(x, "status"))) stop("Cannot record R compiler configuration: ", key)
  paste(x, collapse = "\n")
}
keys <- c("CXX17", "CXX17STD", "CXX17FLAGS", "CPPFLAGS", "LDFLAGS")
flags <- setNames(lapply(keys, config), keys)
# Preserve explicit compiler flags (for example macOS -arch arm64); reject
# wrapper executables whose effective compiler is not identified here.
command <- scan(text = flags$CXX17, what = character(), quiet = TRUE)
if (!length(command) || !grepl("^(clang\\+\\+|((x86_64|aarch64)-linux-gnu-)?g\\+\\+(-[0-9]+)?|c\\+\\+)$", basename(command[[1L]]))) stop("Unsupported compiler wrapper")
compiler <- Sys.which(command[[1L]])
if (!nzchar(compiler)) stop("Cannot locate C++17 compiler")
version <- system2(compiler, "--version", stdout = TRUE, stderr = TRUE)
if (!is.null(attr(version, "status"))) stop("Cannot identify C++17 compiler")
headers <- system.file("include", package = "cpp11", mustWork = TRUE)
paths <- sort(list.files(headers, recursive = TRUE, full.names = TRUE))
header_hashes <- setNames(lapply(paths, function(p) digest::digest(file = p, algo = "sha256")), substring(paths, nchar(headers) + 2L))
user_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = path.expand("~/.R/Makevars"))
site_makevars <- Sys.getenv("R_MAKEVARS_SITE", unset = file.path(R.home("etc"), "Makevars.site"))
# Unrecorded Makevars could silently override the compiler/flags we just read.
for (p in c(user_makevars, site_makevars)) if (file.exists(p) && file.info(p)$size > 0) stop("Native receipt build requires empty user/site Makevars; use explicit reviewed package flags")
jsonlite::write_json(list(schema = 1L, source_revision = unname(read.dcf("DESCRIPTION")[1, "Config/nashtree/SourceRevision"]),
  R = as.character(getRversion()), platform = R.version$platform,
  compiler = normalizePath(compiler), compiler_sha256 = digest::digest(file = compiler, algo = "sha256"),
  compiler_version = version, flags = flags, package_makevars = readLines("src/Makevars"),
  cpp11_version = as.character(packageVersion("cpp11")), cpp11_headers = header_hashes,
  source_inventory_sha256 = digest::digest(file = "inst/native-source-sha256.txt", algo = "sha256")),
  "inst/native-build.json", auto_unbox = TRUE, pretty = TRUE)
