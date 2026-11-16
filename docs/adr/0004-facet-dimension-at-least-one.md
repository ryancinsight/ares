# 0004 — A traction facet needs a spatial dimension of at least one

- **Status:** Accepted
- **Date:** 2026-11-15
- **Driver:** board item `ARES-FACET-DIM-ZERO`, delivered by PR #87
- **Class:** `[major]` — a call that compiled stops compiling, workspace
  0.3.0 under Cargo's leftmost-nonzero rule (the release of ADR 0002).
- **Correction:** the `BREAKING CHANGE` footer of the first commit of PR #87 says
  `TractionFacet` no longer compiles at `D = 0` and that no `D = 0` call ever
  completed work. Both are false; this text and the changelog carry the
  runs below.

## Context

`facet_measure<T, const D: usize>` computes `D - 1` edges. `D` is a
compile-time choice, and nothing bounded it below, so a zero dimension
compiled. At the basis revision, probes of `D = 0` gave:

| Call | Basis | Head |
| ---- | ----- | ---- |
| `facet_measure::<f64, 0>` | panics at the subtraction (debug) | E0080 |
| `try_new` with one facet | panics at the same site | E0080 |
| `try_new(&[], &[[]])` | returns `Ok` | E0080 |
| `TractionFacet::<f64, 0>::new([], [])` | runs | runs |

The release-build panic "at the Gram-matrix index" and the earlier panic of
`add_consistent_loads` in `as_chunks_mut::<0>` come from the board item's
judge run and were not re-run. No run returned a wrong value. The defect was
an undocumented panic, plus a constructor that accepted a boundary it could
not use, and the documentation stated no dimension contract.

## Decision

`D >= 1` is the contract of `facet_measure`, `TractionFacet` and
`TractionBoundary`, stated in their rustdoc. One inline `const` assertion in
`facet_measure` enforces it for the paths that reach that function; for
`TractionFacet` the contract is documentation only, because `new` stores its
arrays and never calls `facet_measure`:

```rust
const { assert!(D >= 1, "a facet needs at least one spatial dimension") }
```

An inline `const` block is evaluated when its function is instantiated, so
the assertion fires for every caller that reaches `facet_measure::<_, 0>`:
`TractionBoundary::try_new`, `add_consistent_loads` and
`ares-coupling`'s `interface_work`. The `try_new` cases were run; the
`add_consistent_loads` and `interface_work` cases follow from the same
instantiation rule and were not run. The error is E0080, raised at
monomorphization: `cargo check` does not report it, `cargo build` and doctests
do. A `compile_fail` doctest on `facet_measure` and one on `TractionBoundary`
pin both entry points.

## Alternatives rejected

**A run-time `assert!` or a typed error.** Rejected: `D` is a type-level
constant, so a run-time check reports at execution a defect the caller fixed
at the call site, and a typed error would make every correct caller handle a
case that cannot occur.

**Assert in each public item.** Rejected: `facet_measure` is the single
computation every measuring path reaches, so one assertion covers every such
caller, including ones added later. The price is that `TractionFacet::new`,
which does no measuring, keeps accepting `D = 0`; a `D = 0` facet cannot be
measured, so it fails to compile at the first path that would use it.

**A `D: NonZero` type parameter.** Rejected: stable Rust has no const-generic
bound, and a marker trait over `const D: usize` would add a seam for one
inequality.

## Consequences

- Breaking: `facet_measure::<_, 0>`, and the `try_new` and
  `add_consistent_loads` instantiations at `D = 0`, no longer compile. They
  panicked at run time, except `try_new(&[], &[[]])`, which returned `Ok`
  and now fails to compile. The migration is `D >= 1`.
- `TractionFacet::<_, 0>::new` still compiles and runs.
- The other `as_chunks::<D>` and `as_chunks_mut::<D>` sites in `ares` and
  `ares-coupling` are outside this decision.
- `cargo-semver-checks` compares rustdoc JSON, whose `Function` item carries
  `sig`, `generics`, `header` and `has_body` but no body (rustdoc-types
  0.61.0, `src/lib.rs:1209-1222`), so the doctest pair, not the SemVer gate,
  settles the break.

## Evidence that would overturn this

A domain that needs a zero-dimensional boundary (a point load on a
one-node system) would justify a distinct type with its own contract, not a
relaxed bound here.
