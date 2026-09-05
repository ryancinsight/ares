# 0001 — The Athena seam is a separate crate

- **Status:** Accepted
- **Date:** 2026-09-04
- **Driver:** atlas ADR 0057 phase A5 (`ATLAS-ARES-PROMOTION-2026-09-03`)
- **Class:** `[arch]` `[patch]` — no published surface exists yet to break.

## Context

Phase A5 requires assembly to an Athena linear system. Athena's seam is

```rust
pub trait LinearOperator<B: KrylovBackend> {
    fn dimension(&self) -> usize;
    fn apply(&self, backend: &B, input: B::View<'_>, output: B::ViewMut<'_>)
        -> Result<(), B::Error>;
}
```

Two properties of that signature decide the structure.

The views are **backend-associated types**. An implementation generic over
`B: KrylovBackend` receives an opaque `B::View<'_>` with no method for reading
element data, so a generic implementation cannot exist — the seam can only be
implemented against a **named** backend.

The error is fixed to **`B::Error`**. An implementation cannot introduce a
failure mode the backend does not already name.

The only host backend is `athena-leto`'s `LetoBackend<T>`, whose views are
`leto::ArrayView1`. `athena-leto` is not `no_std` and links `std` through
`leto`; `athena-core` alone is `no_std` but insufficient, because the views
cannot be read without naming a backend.

`ares` is `#![no_std]`, `#![forbid(unsafe_code)]`, and allocation-free, and
depends only on `aequitas`, `eunomia`, and `proteus` — vocabulary crates with
no infrastructure of their own.

## Decision

Split the repository into a workspace:

- `crates/ares` (`ares-solid`) — the domain core, unchanged in its properties.
- `crates/ares-operator` — the seam, linking `std`, depending on `ares`.

Dependencies run strictly inward. `ares` gains no edge to Athena, Leto, or any
solver.

## Alternatives rejected

**Implement the seam in `ares` and drop `#![no_std]`.** One crate, no
workspace. Rejected: it puts a solver-backend dependency inside the domain
core, which is the coupling ADR 0055's substrate contract exists to prevent,
and it forecloses embedded and kernel-side use of the balance code for a
capability only the host solver needs.

**Gate the seam behind a cargo feature.** Rejected on the stronger ground that
it makes the shipped configuration the one CI does not build by default. A
feature-gated solver path is an untested path, and the failure mode is that the
untested one is the only one users run.

**Return a custom error from `apply`.** Not available: the trait fixes the
error to `B::Error`, which for `LetoBackend` is a closed non-exhaustive enum
with no variant for a degenerate element.

That constraint shaped the design upstream rather than being worked around.
`SimplexMesh::try_new` establishes that every cell integrates — by running the
same `shape_gradients` call assembly will run, not by inferring it from the
measure — and `DirichletConditions::try_new` validates conditions against the
mesh. By the time an operator exists, the only reachable failures are shape
mismatches, which `LetoBackendError::LengthMismatch` names exactly.

## Consequences

The domain core keeps `no_std`, forbids `unsafe`, and allocates nothing. The
seam crate allocates the one scratch buffer the constrained action needs, and
holds it behind a `RefCell` because `apply` takes `&self` while a solver
iterates. That is the only interior mutability in the repository, confined to a
buffer no caller can observe; `ares` keeps the scratch as an explicit parameter
so a caller holding a `&mut` buffer never pays for a borrow flag.

Publication becomes two crates in dependency order, `ares-solid` then
`ares-operator` (atlas ADR 0057 phase A9).

A second backend — Athena's `HephaestusBackend` for accelerators — becomes a
second `impl` in this crate rather than a second crate, since it shares the
domain dependency and the `std` linkage.

## Verification

The architecture test asserts the edge set: `ares` depends on no solver crate,
`ares-operator` depends on `ares`, and there is no edge back. The seam carries
relay-fidelity tests asserting bitwise equality between what Athena receives
and what the domain crate produces, so the adapter cannot silently transform
anything in passing.

## Revisions

**2026-09-04 — the seam crates are named for their concern, not their
dependency.** This record originally named them `ares-athena` and, alongside
it, `ares-harmonia`. Both are the `<host>-<sibling>` shape that AGENTS.md
`standards: Naming prohibition` and `architecture_scoping: Upstream ownership`
prohibit: the name states the dependency rather than the concern, reuses the
sibling's identity for a second thing, and rots when the sibling is renamed or
replaced. They are now `ares-operator` (it presents a linear operator) and
`ares-coupling` (it presents a coupling partition).

The decision this record makes — that the seam is a *separate crate* — is
unchanged, and the reasoning above stands as written. Only the names moved.

Driving evidence: nine crates stack-wide carry the prohibited shape, tracked as
[`atlas#sibling-named-crates`](https://github.com/ryancinsight/atlas/blob/main/backlog.md#sibling-named-crates);
`ares` is the first remediation because both crates were created the same day,
are unpublished, and have no consumer to migrate.

One defect the rename does **not** fix: `ares-coupling` is `publish = true` and
depends on `harmonia`, which is `publish = false`. `engineering_gates: Publish
pipelines` makes that a topology defect, and it is resolved by graduating
`harmonia` rather than in this repository.
