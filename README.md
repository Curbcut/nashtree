# nashtree

nashtree is Curbcut's long-term fork of
[stochtree](https://github.com/StochasticTree/stochtree), maintained for the
NASH+ system of models. It trains BART and BCF models for
[nash](https://github.com/Curbcut/nash) and carries the native forest core
that serves them.

It diverges from upstream in:

- sampler correctness fixes (MCMC proposal probabilities and split extrema, GFR
  NaN ordering, the Gaussian random-effects variance conditional, tree
  reconstruction with recycled node IDs);
- random-intercept support in retention modes and for unseen groups;
- selected and importance-only posterior retention and bulk sample retention;
- batched numeric preprocessing and specialized residual updates;
- the `nashtree_forest_v1` binary format, its R-independent inference core,
  the producer-neutral writer and build receipts (`native/`), and the
  inference-only `nashtreeInference` bindings built from the same core
  (`native/r-package/`).

`NASH-FORK.md` records each change and its evidence.

## Using it

nash pins an exact commit in its `DESCRIPTION` `Remotes` field, and builds its
inference bindings from the same commit, so install nashtree through nash
(`pak::pak("Curbcut/nash")`) rather than on its own. To work on the fork,
install a checkout into an isolated library with `R CMD INSTALL .`.

nashtree and stochtree register the same model classes and S3 methods. Never
load both in one R session; compare them in separate processes.

## Compatibility

The package, DLL and configure identity are `nashtree`. The C++ `StochTree`
namespace, the public R classes and APIs, and historical JSON fields (including
`stochtree_version`) are unchanged, so models saved by the stochtree fork load
in nashtree and serialize identically.

## Upstream

Upstream changes are not merged wholesale. Each one comes in only after a
focused regression review against NASH fits. `upstream-main` mirrors
upstream's `main` for reference. Original author and licence information is
kept in `DESCRIPTION`, `LICENSE` and `inst/COPYRIGHTS`.
