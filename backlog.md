# Ares Backlog

Open items in priority order. The meta item
`ATLAS-ARES-PROMOTION-2026-09-03` (atlas `backlog.md#ares-promotion`) owns
creation, registration and first publication of the repository and is not
copied here. The ladder below follows
[ADR 0003](docs/adr/0003-release-scope.md); every item carries `milestone:`.
Evidence below was read at `basis:` unless stated otherwise.

<a id="M-vessel-wall"></a>
## M-vessel-wall — 0.3.0
- Outcome: a user meshes a thick-walled vessel in Gaia, closes it with Proteus moduli, loads its inner wall through the Harmonia partition, solves through Athena on the Leto host and on the Hephaestus `HostDevice`, and reads displacement, principal stresses and stress invariants; a `Stress`-typed envelope loads a strip through the same partition.
- exit: solve a Gaia-meshed vessel wall to the Lamé solution on host and device — `cargo run --locked -p ares-coupling --example vessel_wall`. Self-checking at `f32` and `f64`: it exits nonzero when the radial displacement, or the radial and hoop stresses read as principal stresses, leave the discretisation bounds derived in its source; when interface work differs from twice the strain energy, or a rigid translation of the configuration changes interface work or stress, beyond the rounding bounds derived there; or when the `HostDevice` solve differs from the host solve beyond the reduction-order bound derived there. `HostDevice` needs no hardware, so no run skips.
- exit: compute principal stresses of a general Cauchy stress — `cargo nextest run --locked --no-tests=fail -p ares-solid --test oracles -E 'test(=constitutive::principal_stresses_of_a_rotated_diagonal_stress_are_its_diagonal)'`; the rotated stress has all six components nonzero.
- exit: compute invariants of a symmetric tensor — `cargo nextest run --locked --no-tests=fail -p ares-solid --test oracles -E 'test(=kinematics::invariants_of_a_rotated_diagonal_tensor_are_its_symmetric_functions)'`.
- exit: solve on the Hephaestus host device to the host answer — `cargo nextest run --locked --no-tests=fail -p ares-operator --test operator -E 'test(=the_host_device_solve_matches_the_host_solve)'`; a `HostDevice` dispatch error fails it.
- exit: load a strip from a typed traction envelope — `cargo nextest run --locked --no-tests=fail -p ares-coupling --test coupling -E 'test(=a_stress_envelope_loads_a_strip_through_the_partition)'`.
- standing: no open `correctness` item; the stack conformance scan's installed-state assertions pass for ares; the ares conformance baseline total stays at or below 1 (atlas `scripts/conformance-baseline.json`).
- needs: none.
- Non-goals: accelerator devices, outbound displacement and stress envelopes, preconditioning, two-way coupling.

<a id="M-typed-exchange"></a>
## M-typed-exchange — 1.0.0
- Outcome: a fluid partition sends `Stress`-typed traction to Ares and receives `Length`-typed displacement and a `Stress`-typed interface stress diagnostic, and the exchange round-trips exactly; the cantilever oracle also runs at `f32`.
- exit: round-trip a traction envelope through export and advance — `cargo nextest run --locked --no-tests=fail -p ares-coupling --test coupling -E 'test(=a_traction_envelope_round_trips_through_export_and_advance)'`.
- exit: exchange typed envelopes with a pressurised strip — `cargo run --locked -p ares-coupling --example typed_exchange`. Self-checking at `f32` and `f64`: the round trip is bitwise exact; the work identity holds within the rounding bound derived in its source; the interface stress envelope matches the strip's closed-form stress, computed from the load and moduli without the solve, within the bound derived there.
- exit: solve a cantilever at `f32` to beam theory — `cargo nextest run --locked --no-tests=fail -p ares-operator --test verification -E 'test(=scalar::the_cantilever_matches_beam_theory_at_f32)'`.
- standing: no open `correctness` item; the ares conformance baseline total is 0; the SemVer gate passes for `ares-solid`, `ares-operator` and `ares-coupling` with `cargo semver-checks --baseline-rev` at the `M-vessel-wall` closing commit its release item names, and against the crates.io release of each published crate; `cargo publish --dry-run --locked` passes for `ares-solid`, then `ares-operator`, against registry providers; README build: every `crates/*/src/lib.rs` carries `#![doc = include_str!("../README.md")]` and `cargo test --doc --locked` passes; `cargo build --examples --locked` and `mdbook test docs/book` pass.
- needs: M-vessel-wall (it extends the typed traction marshalling the strip test introduces).
- Non-goals: two-way coupling, non-conforming transfer, accelerator devices, later-phase physics, publishing `ares-coupling` (it follows `harmonia` graduation).

<a id="ARES-README-CRATE-TABLE"></a>
## ARES-README-CRATE-TABLE — List every workspace crate in the README [patch] — todo
- priority: verification; milestone: M-vessel-wall; needs: none; basis: ab7710580497f09222f3bff6b7518df2e69259e8
- scope: `README.md`
- Defect (read at basis): the `## Crates` table lists `ares-solid` and `ares-operator`; the workspace is `members = ["crates/*"]` and `crates/ares-coupling` exists, which `docs/book/README.md` and `docs/book/stack_position.md` already list as the third crate.
- Outcome: the README table names all three crates with registry names and the ADR 0001 reference.
- Acceptance: every `name` under `crates/*/Cargo.toml` appears in the README table (`cargo metadata --no-deps` names compared with the table rows); `git grep -n 'ares-coupling' README.md` returns the row.
- Next step: add the row, then run the docs and link checks of the committed gate.

<a id="ARES-FACET-DIM-ZERO"></a>
## ARES-FACET-DIM-ZERO — Document and statically reject a zero dimension in facet_measure [major] — todo
- priority: verification; milestone: M-vessel-wall; needs: none; basis: 0ef6b26dc579b47ada1e8a59860ee7a56fd42ff2
- scope: `crates/ares/src/boundary/traction.rs`, `crates/ares/tests/oracles/neumann.rs`, `CHANGELOG.md`, `docs/adr/` (a new ADR for the `D >= 1` contract). Line numbers are basis lines; `ARES-TRACTION-COORDS` rewrote both code files in PR #78 (ADR 0002; at its head `try_new` is :140, `as_chunks_mut` :232), so re-read them at current `main`.
- Defect (run at basis in debug and release by the independent judge of PR 79, not by the filer): public `facet_measure<T, const D: usize>` computes `let edges = D - 1;` (traction.rs:50) with no `D >= 1` guard, and the docs (:5-14, :41-48) state no dimension contract and no `# Panics` section. At `D = 0` a direct call panics at :50 in debug (subtract overflow) and at :68 in release (`edges` wraps, the Gram loop indexes the zero-row `gram`); `TractionBoundary::try_new` with one facet reaches the same sites through :124; `add_consistent_loads` panics first at :183 (`as_chunks_mut::<D>` asserts a non-zero chunk size), so :189 is unreachable. No run returned a wrong value, and `D` is a compile-time choice, so the defect is an undocumented panic (class `verification`). `TractionBoundary::<f64, 0>::try_new(&[], &[[]])` compiles and returns `Ok` at basis.
- Outcome: the rustdoc of `facet_measure`, `TractionFacet` and `TractionBoundary` states the `D >= 1` contract, and one `const { assert!(D >= 1, ...) }` in `facet_measure` rejects `D = 0` at compile time, including through monomorphization of its callers: `TractionBoundary::try_new` (:124), `add_consistent_loads` (:189) and `ares-coupling/src/interface.rs:236` (`interface_work`). Scope is those items only: the other `as_chunks::<D>` / `as_chunks_mut::<D>` sites (`ares/src/assembly/action.rs:66,67,131`, `body.rs:79,80`; `ares-coupling/src/interface.rs:152,178,185,222,223`) are read, not run, and not checked for whether their bounds admit `D = 0`; the next step files each that does.
- Acceptance (facet_measure): a `compile_fail` doctest instantiates `facet_measure::<f64, 0>`; the existing 2-D and 3-D measure oracles keep their assertions.
- Acceptance (constructor): a `compile_fail` doctest calls `TractionBoundary::<f64, 0>::try_new(&[], &[[]])`, which compiles at basis and reaches the assertion through `facet_measure`, so it fails to compile only with the assertion. Stable rustdoc does not check error codes, so each doctest also has a marker-free copy run once by hand that must fail with the const-evaluation error (E0080) while instantiating `facet_measure::<f64, 0>`; record that output in the delivering PR.
- API: public break. The `D = 0` call above compiles at basis and stops compiling, and no documented `D >= 1` contract exists at basis; [major] under the pre-1.0 leftmost-nonzero rule, shipped in the 0.2.0 to 0.3.0 release of `ARES-TRACTION-COORDS`. ADR 0002 records one decision (node borrowing), so this contract gets its own ADR (number claimed in the delivering PR), plus a CHANGELOG migration note: use `D >= 1`; no `D = 0` call ever completed work. `cargo-semver-checks` compares rustdoc, which carries no function bodies (recalled, not looked up), so the doctest pair above is the settling check, not the SemVer gate.
- Prior work: PR #73 (`refs/pull/73/head`) carried `D >= 1` assertions; `ARES-TRACTION-COORDS` left them out of scope.
- Next step: re-read the lines above at current `main`, add the contract text, the assertions and the doctest pair, and file the other `as_chunks` sites that admit `D = 0`.

<a id="ARES-CRATE-README-DOCS"></a>
## ARES-CRATE-README-DOCS — Build each crate README as its crate docs [patch] — todo
- priority: verification; milestone: M-vessel-wall; needs: none; basis: ea52d101fbdaf78fb1f630299983d4d63c100d0c
- scope: `crates/*/README.md`, `crates/*/src/lib.rs`
- Defect (read at basis): `git grep -n include_str -- 'crates/*/src/lib.rs'` has 0 hits, so no gate compiles a README; the install snippets say `0.1.0` (`crates/ares/README.md:19`, `crates/ares-operator/README.md:23`) while the manifest is 0.3.0; `crates/ares/README.md:46` says every oracle runs at `f32` and `f64`, while the Lamé and cantilever oracles (`verification/cylinder.rs`, `structural.rs`) are `f64` only.
- Outcome: each crate's rustdoc front page is its README, the snippets name the manifest version, and the scalar claim matches the tested oracles.
- Acceptance: `git grep -l 'include_str!("../README.md")' -- 'crates/*/src/lib.rs'` lists all 3 crate roots; `cargo test --doc --locked` and `cargo doc --no-deps --locked` pass; `git grep -n '0\.1\.0' -- 'crates/*/README.md'` is empty.
- Next step: merge each crate root's `//!` overview into its README, then replace it with `#![doc = include_str!("../README.md")]`.

<a id="ARES-PACKAGE-CLOSURE"></a>
## ARES-PACKAGE-CLOSURE — Qualify packaged builds against registry providers [patch] — blocked
- priority: verification; milestone: M-vessel-wall; needs: none (external providers); basis: ab7710580497f09222f3bff6b7518df2e69259e8
- scope: `Cargo.toml` provider requirements, `Cargo.lock`, package qualification in CI
- Blocker (crates.io API, 2026-11-13): `proteus-mat` 0.1.1, its latest, requires `eunomia ^0.10.0` while `Cargo.toml` requires eunomia 0.12.0; `leto-ops` `max_stable_version` is 0.46.4 while `Cargo.toml` requires `>=0.47.0, <0.49.0`; `ares-solid`, `ares-operator` and `ares-coupling` are not on the registry. A normalized package build cannot resolve `leto-ops` or the unpublished Ares names from the registry; the eunomia split additionally breaks compilation (PR #73's E0277).
- Re-open trigger: a `proteus-mat` release on eunomia 0.12 and a `leto-ops` 0.47 release are published by their owners. First publication of `ares-solid` is the A9 Ask-User step of atlas `backlog/ares-promotion.md`; `ares-operator` follows it (ADR 0001:82-83); `ares-coupling` waits on `harmonia` graduation (`gap_audit.md#ARES-COUPLING-PUBLISH-CLOSURE`).
- Outcome and acceptance: `cargo package --locked` verification passes for each publishable crate against registry providers, with no provider back-pin and no reduced test gate.
- Next step: when the trigger fires, run `cargo package --locked` for each publishable crate against the registry and record the result here; the E0277 PR #73 reports is not re-run.
