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

<a id="ARES-SPIKE-HOSTDEVICE-CSR"></a>
## ARES-SPIKE-HOSTDEVICE-CSR — Spike: does CsrOperator compose with HostDevice at f32 and f64? [patch] — todo
- priority: feature; milestone: M-vessel-wall; needs: none; basis: 754d6c36f96c1625a93efcd61885c8b9220b0229
- scope: a scratch crate outside the repository (never committed); the deliverable edits `docs/adr/0003-release-scope.md` only.
- Question: athena-hephaestus `CsrOperator` (athena 041afcb `operator.rs:67-72`) bounds `V: DenseVectorOps<D, T>`, `S: SparseOperatorOps<D, T, Matrix = M>` and `T: RealField + Pod + eunomia::Pod`; hephaestus-host implements both traits for `HostDevice` (hephaestus 742e5fb9 `sparse.rs:71`, `dense_vector.rs:124`). Read, not compiled. Does the composition compile and solve at `T = f32` and `T = f64` under the `eunomia` 0.12.0 and `leto` requirements `Cargo.toml` pins (athena Cargo.toml:30 and hephaestus Cargo.toml:68 admit `eunomia >=0.10.0, <0.13.0`; hephaestus Cargo.toml:72 enables `leto` features `parallel` and `mnemosyne-memory`)? Does the lock still resolve, and does it add a second copy of any provider?
- Method: a scratch crate that depends on `athena-hephaestus`, `hephaestus-host` and `hephaestus-core` at those revisions with ares's provider requirements, uploads a small SPD CSR through `CsrOperator::from_parts`, and runs Athena's `Cg` over `HephaestusBackend` on `HostDevice` at both scalars; compiler and resolver output kept verbatim.
- Evidence budget: one scratch crate, at most 3 build attempts (each varying one hypothesis), one solve per scalar. An unresolved lock is a result: record the exact provider version conflict.
- Decision deliverable: rewrite the ADR 0003 "Host device" hypothesis in place with the observed result and the dependency set `ARES-HOSTDEVICE-SOLVE` adds; on a failed composition, apply the first overturning-evidence line of the ADR instead. The ADR stays at or under 150 lines and `adr-index.py check --strict` passes. The rewrite removes this ID from the ADR.
- Next step: add the three Hephaestus crates to the scratch manifest and compile `CsrOperator::from_parts` at `f64`. Links: [ADR 0003](docs/adr/0003-release-scope.md), [M-vessel-wall](#M-vessel-wall).

<a id="ARES-SPIKE-GAIA-ACCESSOR"></a>
## ARES-SPIKE-GAIA-ACCESSOR — Spike: what does a Gaia tetrahedral mesh lend to SimplexMesh? [patch] — todo
- priority: feature; milestone: M-vessel-wall; needs: none; basis: 754d6c36f96c1625a93efcd61885c8b9220b0229
- scope: a scratch crate outside the repository (never committed); the deliverable edits `docs/adr/0003-release-scope.md` only.
- Question: gaia e36d620 (`gaia-mesh` 0.6.10) `TetrahedralMeshBuilder<T>` (`src/domain/mesh/tetrahedral.rs:109`) wraps an `IndexedMesh` and states positive orientation (read, not compiled). Does its output lend `[T; 3]` coordinates and `[usize; 4]` connectivity as slices `ares::SimplexMesh::try_new` accepts, or is a copy needed? Which accessor yields the inner-wall boundary facets a traction needs, how does the example obtain a thick-walled cylinder (a Gaia generator, or `TetrahedralMeshBuilder` calls), and does it build at `f32` and `f64`? Which dependency form does `ares-coupling` use for `gaia-mesh`?
- Method: read the Gaia accessors and rustdoc at that revision, then a scratch probe that builds one tetrahedron at both scalars and passes its output to `SimplexMesh::try_new`; compiler output kept verbatim.
- Evidence budget: one read pass over `src/domain/mesh` and its boundary accessors, one probe crate, at most 3 build attempts.
- Decision deliverable: rewrite the ADR 0003 "Gaia input" hypothesis in place with the accessor shape and the facet source; a needed copy is stated as living in the example, or an upstream Gaia item is filed. The ADR stays at or under 150 lines and `adr-index.py check --strict` passes. The rewrite removes this ID from the ADR.
- Next step: read `tetrahedral.rs` and the `IndexedMesh` accessors, then probe a single tetrahedron. Links: [ADR 0003](docs/adr/0003-release-scope.md), [M-vessel-wall](#M-vessel-wall).

<a id="ARES-SPIKE-STRIP-ENVELOPE"></a>
## ARES-SPIKE-STRIP-ENVELOPE — Spike: how does a strip's traction lay out in a FieldEnvelope? [patch] — todo
- priority: feature; milestone: M-vessel-wall; needs: none; basis: 754d6c36f96c1625a93efcd61885c8b9220b0229
- scope: a scratch crate outside the repository (never committed); the deliverable edits `docs/adr/0003-release-scope.md` only.
- Question: `FieldEnvelope::try_new` takes one `Quantity<T, D>` per cell of a `GridGeometry` and rejects any other count (harmonia 01ebc48 `src/field/envelope.rs:58`); aequitas 75a453a defines `Stress<T>` (`src/systems/si/quantities.rs:81`); `StructuralInterface::read_traction` reads `D` components per facet from a flat slice (`crates/ares-coupling/src/interface.rs:135`). How does a `D`-component traction on the strip fixture's straight interface (`crates/ares-coupling/tests/coupling.rs`, `Fixture`) map onto envelope cells under ADR 0059's node-major order (0059:116-117): one envelope per component, or another shape? Which geometry (shape, spacing, origin, direction) and `Instant` does it carry, and which new `ares-coupling` dependencies follow (`aequitas` moves from dev to normal; `horae` for `Instant`)?
- Method: read `envelope.rs`, `geometry.rs` and ADR 0059 on the node-major order; a scratch probe builds a `Stress`-typed envelope for a strip fixture's facets and reads its values back into the flat order `read_traction` consumes, comparing element by element.
- Evidence budget: one read pass, one probe, at most 3 build attempts.
- Decision deliverable: rewrite the "Typed envelope" hypothesis of ADR 0003 in place with the layout and the signature surface `ARES-STRIP-ENVELOPE-LOAD` implements, or, if no layout carries a strip, apply the ADR's second overturning-evidence line. The ADR stays at or under 150 lines and `adr-index.py check --strict` passes. The rewrite removes this ID from the ADR.
- Next step: read `GridGeometry::try_new` and its direction-matrix validation, then probe a one-column envelope per component. Links: [ADR 0003](docs/adr/0003-release-scope.md), [M-vessel-wall](#M-vessel-wall).

<a id="ARES-SPIKE-SUBSTEP-ADVANCE"></a>
## ARES-SPIKE-SUBSTEP-ADVANCE — Spike: does loading a strip need Partition::advance? [patch] — todo
- priority: feature; milestone: M-vessel-wall; needs: none; basis: 754d6c36f96c1625a93efcd61885c8b9220b0229
- scope: a scratch crate outside the repository (never committed); the deliverable edits `docs/adr/0003-release-scope.md` and, if an upstream constructor is needed, the harmonia board.
- Question: `Substep::new` is `pub(crate)` (harmonia 01ebc48 `src/partition/substep.rs:13`) and `StructuralPartition::advance` takes a `Substep` (`crates/ares-coupling/src/partition.rs:236`), while the inherent `solve_for_traction` takes none. Does `a_stress_envelope_loads_a_strip_through_the_partition` need `advance`, or does `solve_for_traction` suffice? Which public path, if any, produces a `Substep` (Harmonia's coupling driver; harmonia `docs/book/partition_substep.md`)? Its answer mostly serves the `M-typed-exchange` round-trip check, which runs `export` and `advance`.
- Method: read harmonia's driver and the substep book page; a scratch probe drives a stand-in `Partition` through the public driver to a call of `advance`, or fails to construct a `Substep` outside the crate.
- Evidence budget: one read pass, one probe, at most 3 build attempts.
- Decision deliverable: rewrite the Substep sentence of the ADR 0003 "Typed envelope" hypothesis in place; name the entry the strip item uses; if a public constructor is needed, file it on the harmonia board naming the `M-typed-exchange` exit check `a_traction_envelope_round_trips_through_export_and_advance`, which needs it. The ADR stays at or under 150 lines and `adr-index.py check --strict` passes. The rewrite removes this ID from the ADR.
- Next step: read harmonia's driver for every producer of a `Substep`. Links: [ADR 0003](docs/adr/0003-release-scope.md), [M-typed-exchange](#M-typed-exchange).

<a id="ARES-PRINCIPAL-STRESSES"></a>
## ARES-PRINCIPAL-STRESSES — Compute the principal stresses of a Cauchy stress [minor] — todo
- priority: feature; milestone: M-vessel-wall; needs: none; basis: 754d6c36f96c1625a93efcd61885c8b9220b0229
- scope: `crates/ares/src/constitutive/stress.rs` (102 lines; a leaf module under `constitutive/` if the method grows past 500), `crates/ares/tests/oracles/constitutive.rs`, `crates/ares/README.md`, `docs/book/stress.md`, `CHANGELOG.md`.
- Exit check served: compute principal stresses of a general Cauchy stress, `cargo nextest run --locked --no-tests=fail -p ares-solid --test oracles -E 'test(=constitutive::principal_stresses_of_a_rotated_diagonal_stress_are_its_diagonal)'`. At base `git grep -n -i principal_stress 754d6c36f96c1625a93efcd61885c8b9220b0229 -- crates` is empty, so the test is absent and `--no-tests=fail` fails it.
- Outcome: `CauchyStress<T, D>` returns its `D` principal stresses as `Pressure<T>` in a documented order, generic in `T`, for every `D` the crate supports.
- Acceptance: the test rotates a diagonal stress by a rotation that leaves all six components nonzero and asserts the returned values equal the diagonal within a bound derived in the test source (rounding of the rotation product and the eigenvalue sensitivity of a symmetric matrix), at `f32` and `f64` through one generic function; it also asserts the documented order and a repeated-eigenvalue input. No fitted tolerance.
- Next step: derive the eigenvalue method from a source read this session and cite it in the rustdoc, then write the test first. Links: [ADR 0003](docs/adr/0003-release-scope.md), [M-vessel-wall](#M-vessel-wall).

<a id="ARES-STRESS-INVARIANTS"></a>
## ARES-STRESS-INVARIANTS — Compute the invariants of a symmetric tensor [minor] — todo
- priority: feature; milestone: M-vessel-wall; needs: none; basis: 754d6c36f96c1625a93efcd61885c8b9220b0229
- scope: `crates/ares/src/kinematics/tensor.rs` (181 lines; it holds trace, deviator and `double_dot` only), `crates/ares/tests/oracles/kinematics.rs`, `crates/ares/README.md`, `docs/book/kinematics.md`, `CHANGELOG.md`.
- Exit check served: compute invariants of a symmetric tensor, `cargo nextest run --locked --no-tests=fail -p ares-solid --test oracles -E 'test(=kinematics::invariants_of_a_rotated_diagonal_tensor_are_its_symmetric_functions)'`. At base `git grep -n 'fn invariants' 754d6c36f96c1625a93efcd61885c8b9220b0229 -- crates` is empty, so the test is absent and `--no-tests=fail` fails it.
- Outcome: `SymmetricTensor<T, D>` yields `I1`, `I2` and `I3` (ADR 0057:39-40), generic in `T`.
- Acceptance: the test rotates a diagonal tensor with distinct entries and asserts `I1` the sum, `I2` the sum of pairwise products and `I3` the product of those entries, within bounds derived in the test source from the rounding of the rotation product, at `f32` and `f64` through one generic function.
- Gap for the lead: ADR 0003:59 lists the invariants "of a symmetric tensor and of stress", but its mapped exit check covers the tensor only, so a stress form has no oracle and is not in this item; adding it needs a revision note on the ADR and a judge verdict.
- Next step: write the test first, then the three invariants on `SymmetricTensor`. Links: [ADR 0003](docs/adr/0003-release-scope.md), [M-vessel-wall](#M-vessel-wall).

<a id="ARES-HOSTDEVICE-SOLVE"></a>
## ARES-HOSTDEVICE-SOLVE — Solve the constrained stiffness on the Hephaestus HostDevice [minor] — todo
- priority: feature; milestone: M-vessel-wall; needs: ARES-SPIKE-HOSTDEVICE-CSR; basis: 754d6c36f96c1625a93efcd61885c8b9220b0229
- scope: `crates/ares-operator/src/` (a leaf module for CSR assembly and the device solve), `crates/ares-operator/tests/operator.rs`, `crates/ares-operator/Cargo.toml`, `Cargo.toml`, `Cargo.lock`, `docs/adr/0001-athena-seam-as-a-separate-crate.md` (the revision note ADR 0003:43-44 promises in the PR that compiles this route; its :85-87 plans a second `LinearOperator` impl), `docs/book/assembly.md` and `docs/book/solving.md` (both decline an assembled matrix), `crates/ares-operator/README.md`, `CHANGELOG.md`.
- Exit check served: solve on the Hephaestus host device to the host answer, `cargo nextest run --locked --no-tests=fail -p ares-operator --test operator -E 'test(=the_host_device_solve_matches_the_host_solve)'`; a `HostDevice` dispatch error fails it. At base `git grep -n the_host_device_solve 754d6c36f96c1625a93efcd61885c8b9220b0229 -- crates` is empty, so the test is absent and `--no-tests=fail` fails it.
- Outcome: `ares-operator` assembles the constrained stiffness as CSR parts, solves through `CsrOperator` on `HostDevice`, and returns the host solve's displacement at `f32` and `f64`; no run skips.
- Acceptance: the test asserts the device and host displacement fields agree within a bound derived in the test source, through one generic function at both scalars; the host matrix-free solve is the oracle. The first derivation question: two separately terminated Krylov solves may also differ by the convergence tolerance and the operator's conditioning, not only by reduction order.
- Next step: when the spike lands, write the test first, then the CSR assembly. Split the assembly (checked against the matrix-free action) from the device solve if the delivering PR exceeds the 400-line review budget. Links: [ADR 0003](docs/adr/0003-release-scope.md), [M-vessel-wall](#M-vessel-wall).

<a id="ARES-STRIP-ENVELOPE-LOAD"></a>
## ARES-STRIP-ENVELOPE-LOAD — Load a strip from a Stress-typed traction envelope [minor] — todo
- priority: feature; milestone: M-vessel-wall; needs: ARES-SPIKE-STRIP-ENVELOPE, ARES-SPIKE-SUBSTEP-ADVANCE; basis: 754d6c36f96c1625a93efcd61885c8b9220b0229
- scope: `crates/ares-coupling/src/interface.rs`, `crates/ares-coupling/src/partition.rs`, `crates/ares-coupling/tests/coupling.rs`, `crates/ares-coupling/Cargo.toml`, `Cargo.lock`, `crates/ares-coupling/README.md`, `docs/book/coupling.md`, `CHANGELOG.md`.
- Exit check served: load a strip from a typed traction envelope, `cargo nextest run --locked --no-tests=fail -p ares-coupling --test coupling -E 'test(=a_stress_envelope_loads_a_strip_through_the_partition)'`. At base `git grep -n a_stress_envelope_loads 754d6c36f96c1625a93efcd61885c8b9220b0229 -- crates` is empty, so the test is absent and `--no-tests=fail` fails it.
- Outcome: the partition accepts a `Stress`-typed `FieldEnvelope` traction on a strip's interface, in the layout the strip spike fixes, and a traction of another quantity dimension does not compile.
- Acceptance: the test loads the strip fixture of `tests/coupling.rs` from an envelope with a traction that differs per facet and per component, so a wrong cell order fails it, and asserts its displacement equals the displacement from the equivalent flat slice, exactly when the conversion only copies, otherwise within a bound the layout derives; a `compile_fail` doctest passes a quantity of another dimension. The entry (`advance` or `solve_for_traction`) follows the Substep spike.
- Next step: when both spikes land, write the test first, then the envelope-to-slice conversion at the partition boundary. Links: [ADR 0003](docs/adr/0003-release-scope.md), [M-vessel-wall](#M-vessel-wall).

<a id="ARES-VESSEL-WALL-EXAMPLE"></a>
## ARES-VESSEL-WALL-EXAMPLE — Solve a Gaia-meshed vessel wall to the Lame solution on host and device [minor] — todo
- priority: feature; milestone: M-vessel-wall; needs: ARES-SPIKE-GAIA-ACCESSOR, ARES-PRINCIPAL-STRESSES, ARES-HOSTDEVICE-SOLVE; basis: 754d6c36f96c1625a93efcd61885c8b9220b0229
- scope: `crates/ares-coupling/examples/vessel_wall.rs`, `crates/ares-coupling/Cargo.toml`, `Cargo.toml`, `Cargo.lock`, `README.md`, a book chapter that includes the example source (`docs/book/SUMMARY.md`), `CHANGELOG.md`.
- Exit check served: solve a Gaia-meshed vessel wall to the Lame solution on host and device, `cargo run --locked -p ares-coupling --example vessel_wall`. At base `git ls-tree -r --name-only 754d6c36f96c1625a93efcd61885c8b9220b0229 | grep -c examples/` prints 0, so the command fails.
- Outcome: the example meshes a thick-walled vessel in Gaia, loads its inner wall through the partition, solves on the host and on `HostDevice`, and reads displacement and principal stresses.
- Acceptance: it exits nonzero, at `f32` and `f64`, when the radial displacement or the radial and hoop stresses read as principal stresses leave bounds derived in its source; when interface work differs from twice the strain energy, or a rigid translation changes interface work or stress, beyond derived rounding bounds; or when the `HostDevice` solve differs from the host solve beyond the derived reduction-order bound. The `f32` stress bound is derived in the source from the scalar's epsilon and the mesh error, with no fitted single-mesh tolerance (`crates/ares-operator/tests/verification/cylinder.rs:25-32`, `scalar.rs:18-27`).
- Gap for the delivering PR, which neither the ladder nor ADR 0003 fixes: the vessel dimensions and loading, the axial constraint that makes a plane-strain closed form apply to a 3-D Gaia mesh, the refinement the `f32` bound uses, and how the `HostDevice` solve receives the partition-built load (`solve_for_traction` is fixed to `Cg::<LetoBackend<T>>` and the partition exposes no load). The existing closed form is 2-D plane strain (`verification/cylinder.rs`); the example header records each choice.
- Next step: when the spike and both items land, record those choices, write the host check first, then the device and `f32` checks. Split at the 400-line review budget, the Gaia adapter and host Lame check first. Links: [ADR 0003](docs/adr/0003-release-scope.md), [M-vessel-wall](#M-vessel-wall).

<a id="ARES-PACKAGE-CLOSURE"></a>
## ARES-PACKAGE-CLOSURE — Qualify packaged builds against registry providers [patch] — blocked
- priority: verification; milestone: M-vessel-wall; needs: none (external providers); basis: ab7710580497f09222f3bff6b7518df2e69259e8
- scope: `Cargo.toml` provider requirements, `Cargo.lock`, package qualification in CI
- Blocker (crates.io API, 2026-11-13): `proteus-mat` 0.1.1, its latest, requires `eunomia ^0.10.0` while `Cargo.toml` requires eunomia 0.12.0; `leto-ops` `max_stable_version` is 0.46.4 while `Cargo.toml` requires `>=0.47.0, <0.49.0`; `ares-solid`, `ares-operator` and `ares-coupling` are not on the registry. A normalized package build cannot resolve `leto-ops` or the unpublished Ares names from the registry; the eunomia split additionally breaks compilation (PR #73's E0277).
- Re-open trigger: a `proteus-mat` release on eunomia 0.12 and a `leto-ops` 0.47 release are published by their owners. First publication of `ares-solid` is the A9 Ask-User step of atlas `backlog/ares-promotion.md`; `ares-operator` follows it (ADR 0001:82-83); `ares-coupling` waits on `harmonia` graduation (`gap_audit.md#ARES-COUPLING-PUBLISH-CLOSURE`).
- Outcome and acceptance: `cargo package --locked` verification passes for each publishable crate against registry providers, with no provider back-pin and no reduced test gate.
- Next step: when the trigger fires, run `cargo package --locked` for each publishable crate against the registry and record the result here; the E0277 PR #73 reports is not re-run.
