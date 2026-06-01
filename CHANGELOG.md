# Changelog

## Unreleased

### Fixed
- `xy_ring(xx::Vector{String})` no longer errors on an undefined `kwargs` reference.
- `is_rkmtx` now iterates the interior grid correctly (was `for (i,j) in (1:n-1,1:m-1)`,
  which never visited the `(i,j)` pairs); added a regression test.
- `xy_ring(n, m; varnames=...)` now honors `varnames` instead of hardcoding `x`/`y`.

### Changed
- Coefficient-ring keyword renamed `coef_ring` → `coeff` (consistent with
  SemistandardTableaux).
- `extract_vars` returns a concrete `elem_type(R)[]` vector (was abstract
  `MPolyRingElem[]`).
- `partition2drift` gains a keyword form `partition2drift(lambda; ff=, n=)`; the
  positional `partition2drift(lambda, ff[, n])` form is retained for back-compatibility.
- `all_bpds` is now re-exported (matching the existing `BPD`/`draw_bpd` re-exports).

### Added
- Docstrings for `drift2bin`, `drift_poly`, `partition2drift`, `bpd2drift`, `markconfig`.
- Integer drift-entry encoding legend documented at the `Drift` type definition.
- `test/drift_polys.jl` (re-enabled in the test list; fixed `draw_drifts.jl` →
  `draw_drift.jl` filename typo).
- `[extras]`/`[targets]` declaring `Test` so `Pkg.test()` resolves on Julia 1.12.

### Notes
- Canonical `Drift` field is `.mtx`.
- Drift-marking pipeline kept at **1d-min**: `dc2sd`/`schub_drifts`/`tableau_components`
  remain commented out, pending decoupling from `schur_poly` (1d-full follow-up).
