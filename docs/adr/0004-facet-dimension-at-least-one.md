# 0004 — A traction facet needs a spatial dimension of at least one

- **Status:** Accepted
- **Date:** 2026-11-15
- **Driver:** board item `ARES-FACET-DIM-ZERO`, delivered by PR #87
- **Class:** `[major]` — a call that compiled stops compiling, workspace
  0.3.0 under Cargo's leftmost-nonzero rule (the release of ADR 0002).

## Context

`facet_measure<T, const D: usize>` computes `D - 1` edges. `D` is a
compile-time choice, and nothing bounded it below, so a zero dimension
compiled. Run at the basis revision, `D = 0` panicked at the subtraction in a
debug build and at the Gram-matrix index in a release build; reached through
`TractionBoundary::try_new` or `add_consistent_loads` it panicked at the same
sites or, for the latter, earlier in `as_chunks_mut::<0>`. No run returned a
wrong value. The defect was an undocumented panic, and the documentation
stated no dimension contract.

## Decision

`D >= 1` is the contract of `facet_measure`, `TractionFacet` and
`TractionBoundary`, stated in their rustdoc. One inline `const` assertion in
`facet_measure` enforces it:

```rust
const { assert!(D >= 1, "a facet needs at least one spatial dimension") }
```

An inline `const` block is evaluated when its function is instantiated, so
the assertion fires for every caller that reaches `facet_measure::<_, 0>`:
`TractionBoundary::try_new`, `add_consistent_loads` and
`ares-coupling`'s `interface_work`. The error is E0080, raised at
monomorphization: `cargo check` does not report it, `cargo build` and doctests
do. A `compile_fail` doctest on `facet_measure` and one on `TractionBoundary`
pin both entry points.

## Alternatives rejected

**A run-time `assert!` or a typed error.** Rejected: `D` is a type-level
constant, so a run-time check reports at execution a defect the caller fixed
at the call site, and a typed error would make every correct caller handle a
case that cannot occur.

**Assert in each public item.** Rejected: `facet_measure` is the single
computation every dimension-dependent path reaches, so one assertion covers
every caller, including ones added later.

**A `D: NonZero` type parameter.** Rejected: stable Rust has no const-generic
bound, and a marker trait over `const D: usize` would add a seam for one
inequality.

## Consequences

- Breaking: a `D = 0` instantiation, which compiled and could not complete
  work, no longer compiles. The migration is `D >= 1`.
- The other `as_chunks::<D>` and `as_chunks_mut::<D>` sites in `ares` and
  `ares-coupling` are outside this decision.
- `cargo-semver-checks` compares rustdoc, which carries no function bodies,
  so the doctest pair, not the SemVer gate, settles the break.

## Evidence that would overturn this

A domain that needs a zero-dimensional boundary (a point load on a
one-node system) would justify a distinct type with its own contract, not a
relaxed bound here.
