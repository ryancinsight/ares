# 0004 — Traction, mesh and interface types need a spatial dimension of at least one

- **Status:** Accepted
- **Date:** 2026-11-15
- **Driver:** board item `ARES-FACET-DIM-ZERO`, delivered by PR #87
- **Class:** `[major]` — a call that compiled stops compiling, workspace
  0.3.0 under Cargo's leftmost-nonzero rule (the release of ADR 0002).
- **Revision:** 2026-11-24, PR #102 (`ARES-ASSEMBLY-DIM-ZERO`): the contract
  extends to `SimplexMesh` and `StructuralInterface`; the facet decision is
  unchanged.
- **Revision:** 2026-11-27, 5643cd0ea13d (PR #106): the title named only a
  traction facet; it now reads "Traction, mesh and interface types need a
  spatial dimension of at least one", which spans `TractionFacet`,
  `TractionBoundary`, `SimplexMesh` and `StructuralInterface`, the types the
  contract has covered since PR #102. The decision is unchanged.
- **Correction:** the `BREAKING CHANGE` footer of the commit titled
  `fix(ares): Reject a zero facet dimension` says `TractionFacet` no longer
  compiles at `D = 0` and that no `D = 0` call ever completed work. Both are
  false; this text and the changelog carry the runs below.

## Context

`facet_measure<T, const D: usize>` computes `D - 1` edges. `D` is a
compile-time choice, and nothing bounded it below, so a zero dimension
compiled. Probes of `D = 0`, run by the PR's judge, gave:

| Call | Before PR #87 | After PR #87 |
| ---- | ------------- | ------------ |
| `facet_measure::<f64, 0>` | panics at the subtraction (debug) | E0080 |
| `try_new` with one facet | panics at the same site | E0080 |
| `try_new(&[], &[[]])` | returns `Ok` | E0080 |
| `TractionFacet::<f64, 0>::new([], [])` | runs | runs |

The release-build panic "at the Gram-matrix index" and the earlier panic of
`add_consistent_loads` in `as_chunks_mut::<0>` come from the board item's
judge run and were not re-run. No run returned a wrong value. The defect was
an undocumented panic, a constructor that returned `Ok` for
`try_new(&[], &[[]])` at `D = 0`, and no stated dimension contract.

## Decision

`D >= 1` is the contract of `facet_measure`, `TractionFacet` and
`TractionBoundary`, stated in their rustdoc. One inline `const` assertion in
`facet_measure` enforces it for the paths that reach that function; for
`TractionFacet` the contract is documentation only, because `new`, the
accessors and the derived traits never call `facet_measure`:

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
do. A `compile_fail` doctest on `facet_measure` pins the direct call, and one
on `TractionBoundary` pins the `try_new` path.

## Alternatives rejected

**A run-time `assert!` or a typed error.** Rejected: `D` is a type-level
constant, so a run-time check reports at execution a defect the caller fixed
at the call site, and a typed error would make every correct caller handle a
case that cannot occur.

**Assert in each public item.** Rejected: `facet_measure` is the single
computation every measuring path reaches, so one assertion covers every such
caller, including ones added later. The price is that `TractionFacet::new`,
its accessors and its derived traits never reach `facet_measure` and keep
accepting `D = 0`. Only the paths that reach `facet_measure`, namely
`TractionBoundary::try_new`, `add_consistent_loads` and `interface_work`, fail
to compile; a `D = 0` facet can still be built, read, copied and compared.

**A `D: NonZero` type parameter.** Rejected: stable Rust has no const-generic
bound, and a marker trait over `const D: usize` would add a seam for one
inequality.

## Consequences

- Breaking: `facet_measure::<_, 0>`, and the `try_new` and
  `add_consistent_loads` instantiations at `D = 0`, no longer compile. They
  panicked at run time, except `try_new(&[], &[[]])`, which returned `Ok`
  and now fails to compile. The migration is `D >= 1`.
- `TractionFacet::<_, 0>::new` still compiles and runs.
- The `as_chunks::<D>` and `as_chunks_mut::<D>` sites of assembly and coupling
  are covered by the extension below.
- `cargo-semver-checks` compares rustdoc JSON, whose `Function` item carries
  `sig`, `generics`, `header`, `has_body` and `default_unstable` but no body
  (rustdoc-types 0.61.0, `src/lib.rs:1209-1224`), so the doctest pair, not the
  SemVer gate, settles the break.

## Evidence that would overturn this

A domain that needs a zero-dimensional boundary (a point load on a
one-node system) would justify a distinct type with its own contract, not a
relaxed bound here.

## Extension to assembly and coupling

`D >= 1` is also the contract of `SimplexMesh` and `StructuralInterface`. A
field carries `D` components per node and is read as `D`-wide chunks, and a
zero width has no nodal view. Each constructor holds one inline `const`
assertion of the same form as `facet_measure`, with its own message, so the
E0080 names the type: `SimplexMesh::try_new` ("a mesh needs at least one
spatial dimension") and `StructuralInterface::try_new` ("an interface needs at
least one spatial dimension").

At `D = 0` the calls ran as follows. Both constructors returned `Ok`, so the
failure surfaced only at the first operation. The five panics ran in one debug
build, each under `catch_unwind`, with the message `chunk size must be
non-zero`; `interface_work` ran in a separate build.

| Call at `D = 0` | Before PR #102 | After PR #102 |
| --- | --- | --- |
| `SimplexMesh::try_new(&[[]], &[[0]])` | `Ok` | E0080 |
| `StructuralInterface::try_new(&[0], &[[]], &mesh)` | `Ok` | E0080 |
| `add_body_force` | panics, `body.rs:79` | E0080 at `try_new` |
| `internal_forces` | panics, `action.rs:66` | E0080 at `try_new` |
| `cell_strains` | panics, `action.rs:131` | E0080 at `try_new` |
| `read_traction` | panics, `interface.rs:153` | E0080 at `try_new` |
| `write_displacement` | panics, `interface.rs:179` | E0080 at `try_new` |
| `interface_work` | E0080 in `facet_measure`, `interface.rs:237` | E0080 at `try_new` |

Every safe path to the operations above goes through one of the two
constructors, because both types keep their fields private. The adjacent
`D = 1` bar is retained as value oracles, `oracles::bar` for the mesh and
`a_point_interface_on_a_bar_exchanges_exactly` for the interface: a unit bar
under a body force of 2 assembles the exact nodal loads `[1, 1]`.

The interface assertion is redundant for a caller that obtains its mesh from
`SimplexMesh::try_new`, because that constructor already rejects `D = 0`. It
is kept for the diagnostic and so the invariant holds if the mesh bound is
relaxed; its `compile_fail` example builds no mesh, so only the interface
assertion can fail it.

**Alternatives rejected.** A typed error: `D` is a type-level constant, so a
check at run time reports a defect the caller fixed at the call site. A check
at each split: it defers a condition already known from the const argument.
Point assembly: it needs the field and interface semantics defined for every
operation above, and this decision defines none.

**Consequences.** Breaking, `[major]` before 1.0:
`SimplexMesh::<_, 0, N>::try_new` and `StructuralInterface::<0>::try_new` no
longer compile. The migration is `D >= 1`. The E0080 is raised at
instantiation, as for `facet_measure` above, so the two `compile_fail`
doctests, not `cargo check` or the SemVer gate, pin the break.
