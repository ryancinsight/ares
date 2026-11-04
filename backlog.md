# Ares backlog

The active target is 0.3.0. Acceptance checklists live with their canonical
items; no separate root checklist is introduced.

<a id="ares-traction-domain-20261104"></a>
## ARES-TRACTION-DOMAIN-20261104 — Retain validated coordinates [major] — blocked

- Priority: correctness; owner: codex:01a17995-4c0b-701d-8625-8650165230c2; target: 0.3.0.
- Basis: construction validates one node slice but the old assembly API accepts
  another. Both f32/f64 native witnesses show a shorter-array panic and silently
  zero load from a collapsed replacement. This item owns traction identity,
  not the facet-measure numerical algorithm or another repository's migration.
- [x] Reproduce both public behaviors before implementation and retain the
  unchanged-source receipts; strengthen mutation/lifetime compile contracts.
- [x] Review the immutable-borrow design, record
  [ADR 0002](docs/adr/0002-traction-coordinate-domain.md), and migrate all local
  operator, coupling, oracle and book callers without a compatibility shim.
- [x] Pass the focused generic 1-D/2-D/3-D, additive, shape-error and geometry
  reconstruction tests at f32/f64/F32/F64 and the five new compile-fail cases.
- [ ] Complete final-revision full/minimal/release/MSRV/API/package/book/cost
  checks, independent source review, normal hooks and hosted qualification.
- [ ] Resolve the release prerequisites below before merge/publication; verify
  published source and permitted documentation deployment after delivery.

<a id="ares-package-closure-20261104"></a>
## ARES-PACKAGE-CLOSURE-20261104 — Coherent published providers [major] — blocked

- Priority: correctness; owner: existing provider/cohort owners with Ares
  coordination by codex:01a17995-4c0b-701d-8625-8650165230c2; target: 0.3.0.
- Basis: the normalized baseline core archive builds, but its tests fail with
  E0277 in `tests/oracles/support.rs`: `T: eunomia::RealField` from 0.12 does
  not imply Proteus 0.1.1's Eunomia 0.10 trait. This is reproduced on unchanged
  source and is not introduced by traction identity.
- [x] Preserve the failed baseline package-test receipt separately from source
  correctness tests and identify the existing [Proteus #94](https://github.com/ryancinsight/proteus/pull/94)
  provider/cohort transition rather than adopt its unmerged branch.
- [ ] Existing owners qualify and deliver compatible latest provider and
  consumer requirements/locks; preserve Ares #69 dependency ownership.
- [ ] Qualify every normalized Ares workspace component against that published
  closure before release. No provider back-pin or reduced test gate is allowed.
- First-publication state is separate: all three registry names returned 404
  on 2026-11-04. The pinned pending-version publisher skips never-published
  names; no authentication, token, publisher or permission change is attempted.
