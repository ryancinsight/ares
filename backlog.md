# Ares Backlog

Open items in priority order. The meta item
`ATLAS-ARES-PROMOTION-2026-09-03` (atlas `backlog.md#ares-promotion`) owns
creation, registration and first publication of the repository and is not
copied here. Items carry no `milestone:` until `ARES-RELEASE-LADDER` lands a
ladder. Evidence below was read at `basis:` unless stated otherwise.

<a id="ARES-TRACTION-COORDS"></a>
## ARES-TRACTION-COORDS — Bind traction loads to the validated coordinates [major] — todo
- priority: correctness; needs: none; basis: ab7710580497f09222f3bff6b7518df2e69259e8
- scope: `crates/ares/src/boundary/traction.rs`, its callers `crates/ares-coupling/src/partition.rs`, `crates/ares-operator/tests/`, `crates/ares/tests/oracles/neumann.rs`, `docs/book/boundary.md`, `CHANGELOG.md`, `docs/adr/`
- Defect (confirmed by reading `traction.rs` at basis): `TractionBoundary::try_new` validates facets against one `nodes` slice and retains only `facets`; `add_consistent_loads(&self, nodes, loads)` takes a second slice. Its shape check uses the replacement's `nodes.len() * D`, then `nodes[*node]` panics when the replacement is shorter than a facet index, and coincident replacement nodes give facet measure 0, return `Ok` and add no load.
- Outcome: the boundary retains the validated coordinates by immutable borrow and `add_consistent_loads` takes only `loads`; every in-repo caller migrates in the same change, no forwarding method.
- Acceptance: native tests at f32 and f64 (one-node replacement; coincident replacement) that fail at the basis and pass with the binding; existing force, moment and energy oracles keep their assertions.
- API: public break, workspace 0.2.0 to 0.3.0 (pre-1.0 leftmost-nonzero rule), so it needs ADR 0002 (absent at basis, `docs/adr` holds only 0001) and a SemVer-gate [major] classification with migration note.
- Prior work: released PR #73 (`fix/traction-mesh-identity`, kept at `refs/pull/73/head`, label `released`) carries the fix, ADR 0002 and `tests/traction_identity.rs`. Its recorded claims (witnesses fail on unchanged source; 108 tests pass in hosted CI) are not re-run here.
- Next step: `git fetch origin pull/73/head:fix/traction-mesh-identity`, port onto current main without the PR's own `backlog.md`/`gap_audit.md` (this board replaces them), confirm the witnesses fail at the basis, then run the committed gate. Publication stays held by `ARES-PACKAGE-CLOSURE`.

<a id="ARES-RELEASE-LADDER"></a>
## ARES-RELEASE-LADDER — Draft the release scope and milestone ladder [patch] — todo
- priority: verification; needs: none; basis: ab7710580497f09222f3bff6b7518df2e69259e8
- scope: `docs/adr/` (release-scope ADR; ADR 0002 is reserved by `ARES-TRACTION-COORDS`), `docs/adr/README.md`, `backlog.md`
- Outcome: release-scope ADR and milestone ladder per the stack roadmap rule (sprint: roadmap). The ladder opens this board; every item then carries `milestone:`.
- Acceptance: the ADR records each capability area of `ares-solid`, `ares-operator` and `ares-coupling` as delivered, in scope or a non-goal; each milestone exit check is a runnable check that fails without its feature; an independent judge verdict maps every capability area and open feature item to a check or non-goal.
- Inputs: `README.md` scope paragraph, atlas ADRs 0055, 0057 and 0059, `docs/book`. Phase 0 non-goals (plasticity, contact, finite deformation, dynamics, fracture) come from the meta item.
- Next step: draft the ADR from those inputs with a recommended scope, name any scope no source determines as a veto point, and file the ladder under its own `judge-required` PR.

<a id="ARES-README-CRATE-TABLE"></a>
## ARES-README-CRATE-TABLE — List every workspace crate in the README [patch] — todo
- priority: verification; needs: none; basis: ab7710580497f09222f3bff6b7518df2e69259e8
- scope: `README.md`
- Defect (read at basis): the `## Crates` table lists `ares-solid` and `ares-operator`; the workspace is `members = ["crates/*"]` and `crates/ares-coupling` exists, which `docs/book/README.md` and `docs/book/stack_position.md` already list as the third crate.
- Outcome: the README table names all three crates with registry names and the ADR 0001 reference.
- Acceptance: every `name` under `crates/*/Cargo.toml` appears in the README table (`cargo metadata --no-deps` names compared with the table rows); `git grep -n 'ares-coupling' README.md` returns the row.
- Next step: add the row, then run the docs and link checks of the committed gate.

<a id="ARES-FACET-DIM-ZERO"></a>
## ARES-FACET-DIM-ZERO — Document and statically reject a zero dimension in facet_measure [major] — todo
- priority: verification; needs: ARES-TRACTION-COORDS; basis: 0ef6b26dc579b47ada1e8a59860ee7a56fd42ff2
- scope: `crates/ares/src/boundary/traction.rs`, `crates/ares/tests/oracles/neumann.rs`, `CHANGELOG.md`, `docs/adr/` (a new ADR for the `D >= 1` contract). Line numbers are basis lines; `ARES-TRACTION-COORDS` rewrites both code files (at its head `try_new` is :140, `as_chunks_mut` :232), so re-read them once it lands.
- Defect (run at basis in debug and release by the independent judge of PR 79, not by the filer): public `facet_measure<T, const D: usize>` computes `let edges = D - 1;` (traction.rs:50) with no `D >= 1` guard, and the docs (:5-14, :41-48) state no dimension contract and no `# Panics` section. At `D = 0` a direct call panics at :50 in debug (subtract overflow) and at :68 in release (`edges` wraps, the Gram loop indexes the zero-row `gram`); `TractionBoundary::try_new` with one facet reaches the same sites through :124; `add_consistent_loads` panics first at :183 (`as_chunks_mut::<D>` asserts a non-zero chunk size), so :189 is unreachable. No run returned a wrong value, and `D` is a compile-time choice, so the defect is an undocumented panic (class `verification`). `TractionBoundary::<f64, 0>::try_new(&[], &[[]])` compiles and returns `Ok` at basis.
- Outcome: the rustdoc of `facet_measure`, `TractionFacet` and `TractionBoundary` states the `D >= 1` contract, and one `const { assert!(D >= 1, ...) }` in `facet_measure` rejects `D = 0` at compile time, including through monomorphization of its callers: `TractionBoundary::try_new` (:124), `add_consistent_loads` (:189) and `ares-coupling/src/interface.rs:236` (`interface_work`). Scope is those items only: the other `as_chunks::<D>` / `as_chunks_mut::<D>` sites (`ares/src/assembly/action.rs:66,67,131`, `body.rs:79,80`; `ares-coupling/src/interface.rs:152,178,185,222,223`) are read, not run, and not checked for whether their bounds admit `D = 0`; the next step files each that does.
- Acceptance (facet_measure): a `compile_fail` doctest instantiates `facet_measure::<f64, 0>`; the existing 2-D and 3-D measure oracles keep their assertions.
- Acceptance (constructor): a `compile_fail` doctest calls `TractionBoundary::<f64, 0>::try_new(&[], &[[]])`, which compiles at basis and reaches the assertion through `facet_measure`, so it fails to compile only with the assertion. Stable rustdoc does not check error codes, so each doctest also has a marker-free copy run once by hand that must fail with the const-evaluation error (E0080) while instantiating `facet_measure::<f64, 0>`; record that output in the delivering PR.
- API: public break. The `D = 0` call above compiles at basis and stops compiling, and no documented `D >= 1` contract exists at basis; [major] under the pre-1.0 leftmost-nonzero rule, shipped in the 0.2.0 to 0.3.0 release of `ARES-TRACTION-COORDS`. ADR 0002 records one decision (node borrowing), so this contract gets its own ADR (number claimed in the delivering PR), plus a CHANGELOG migration note: use `D >= 1`; no `D = 0` call ever completed work. `cargo-semver-checks` compares rustdoc, which carries no function bodies (recalled, not looked up), so the doctest pair above is the settling check, not the SemVer gate.
- Prior work: PR #73 (`refs/pull/73/head`) carried `D >= 1` assertions; `ARES-TRACTION-COORDS` left them out of scope.
- Next step: after `ARES-TRACTION-COORDS` lands, re-read the lines above, add the contract text, the assertions and the doctest pair, and file the other `as_chunks` sites that admit `D = 0`.

<a id="ARES-PACKAGE-CLOSURE"></a>
## ARES-PACKAGE-CLOSURE — Qualify packaged builds against registry providers [patch] — blocked
- priority: verification; needs: none (external providers); basis: ab7710580497f09222f3bff6b7518df2e69259e8
- scope: `Cargo.toml` provider requirements, `Cargo.lock`, package qualification in CI
- Blocker (crates.io API, 2026-11-13): `proteus-mat` 0.1.1, its latest, requires `eunomia ^0.10.0` while `Cargo.toml` requires eunomia 0.12.0; `leto-ops` `max_stable_version` is 0.46.4 while `Cargo.toml` requires `>=0.47.0, <0.49.0`; `ares-solid`, `ares-operator` and `ares-coupling` are not on the registry. A normalized package build cannot resolve `leto-ops` or the unpublished Ares names from the registry; the eunomia split additionally breaks compilation (PR #73's E0277).
- Re-open trigger: a `proteus-mat` release on eunomia 0.12 and a `leto-ops` 0.47 release are published by their owners. First publication of `ares-solid` is the A9 Ask-User step of atlas `backlog/ares-promotion.md`; `ares-operator` and `ares-coupling` are outside A9 while `harmonia` is `publish = false`.
- Outcome and acceptance: `cargo package --locked` verification passes for each publishable crate against registry providers, with no provider back-pin and no reduced test gate.
- Next step: when the trigger fires, run `cargo package --locked` for each publishable crate against the registry and record the result here; the E0277 PR #73 reports is not re-run.
