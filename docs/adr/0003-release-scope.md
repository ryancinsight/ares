# 0003 — Release scope for 1.0

- **Status:** Proposed
- **Date:** 2026-11-14
- **Driver:** board item `ARES-RELEASE-LADDER`, delivered by PR #85
- **Class:** `[patch]` — records scope and the milestone ladder; changes no
  code.

## Context

Ares has no release tag; the manifest is at 0.3.0. Sources: atlas ADRs 0055,
0057 (Accepted) and 0059 (Proposed), ADRs 0001 and 0002, the READMEs, the book
and the board. ADR 0057 fixes the solve as "backend-neutral over Leto and
Hephaestus" (0057:53-54), the exchange as "typed physical fields per ADR 0050"
(0057:85-88), and every oracle at `f32` and `f64` (condition 5). ADR
0001:82-83 fixes publication. No source states what 1.0 is.

## Decision

1.0 is the complete Phase 0 vertical slice of ADR 0057 plus the one-way
coupling of ADR 0059. A user meshes a solid in Gaia, closes it with Proteus
moduli, loads it directly or through a Harmonia partition with `Stress`-typed
traction envelopes, solves through Athena on the Leto host or a Hephaestus
device, and reads displacement, principal stresses and stress invariants.
Later phases are not 1.0 scope.

The ladder in `backlog.md` has two milestones: `M-vessel-wall` (0.3.0) and
`M-typed-exchange` (1.0.0). With no release tag, the first block takes the
manifest version. The first milestone runs one thin real path through every
seam the scope names: Gaia mesh, Proteus closure, Harmonia partition with a
typed traction envelope, Athena on the Leto host, and Athena on the Hephaestus
host reference device. Device feasibility, the largest open risk, is retired
first. The second milestone completes the typed exchange: displacement and
stress envelopes out, and the round trip through `export` and `advance`. It
comes last because it may need an upstream Harmonia constructor (Hypotheses).

**Hephaestus route.** Ares supplies the constrained stiffness as CSR parts to
athena-hephaestus `CsrOperator`, the only Athena operator over
`HephaestusBackend` (athena 041afcb `operator.rs:67`); the host matrix-free
solve is the oracle. ADR 0001:85-87 planned a second `LinearOperator` impl in
`ares-operator`; `CsrOperator` is generic over every `SparseOperatorOps`
device, so none is needed. CSR assembly allocates and lives in
`ares-operator`, so ADR 0001's "rather than a second crate" holds; ADR 0001
gets its revision note in the PR that compiles this route.

**Publication** follows ADR 0001:82-83: `ares-solid`, then `ares-operator`.
`ares-coupling` publishes once `harmonia` graduates, which the ADR 0001
revision of 2026-09-04 places outside this repository
(`gap_audit.md#ARES-COUPLING-PUBLISH-CLOSURE`). crates.io API, 2026-11-14:
`athena-core` and `athena-leto` 0.3.3, `leto` 0.47.0, `leto-ops` 0.46.4
against a required `>=0.47`; `ares-solid` and `ares-operator` absent. Both dry
runs therefore wait on the providers `ARES-PACKAGE-CLOSURE` names.

## Capability areas

| Crate | Area | State |
| --- | --- | --- |
| ares-solid | Displacement gradient, small strain, symmetric tensor algebra | delivered (`oracles::kinematics`) |
| ares-solid | Invariants `I1`, `I2`, `I3` of a symmetric tensor and of stress (ADR 0057:39-40); `tensor.rs` has trace, deviator and `double_dot` only | in scope, `M-vessel-wall`: test `kinematics::invariants_of_a_rotated_diagonal_tensor_are_its_symmetric_functions` |
| ares-solid | Cauchy stress, mean and deviatoric stress, typed `J2`, von Mises | delivered (`oracles::constitutive`) |
| ares-solid | Principal stresses | in scope, `M-vessel-wall`: test `constitutive::principal_stresses_of_a_rotated_diagonal_stress_are_its_diagonal`; `vessel_wall` reads the Lamé radial and hoop stresses through it |
| ares-solid | Isotropic Hooke over `proteus` moduli; no material constants | delivered (`constitutive::the_crate_names_no_material`) |
| ares-solid | Equilibrium residual, matrix-free internal forces | delivered (`oracles::patch_test`, `oracles::assembled_operator`) |
| ares-solid | Linear simplices in 2-D and 3-D, exact one-point integration | delivered (`oracles::element_*`); for linear simplices, ADR 0057's isoparametric mapping and Gauss quadrature reduce to this |
| ares-solid | Higher-order elements, general quadrature rules | non-goal |
| ares-solid | Mesh validation over a borrowed `SimplexMesh` | delivered (`oracles::mesh_validation`) |
| ares-solid | Gaia tetrahedral mesh as solve input | in scope, `M-vessel-wall`: example `vessel_wall` |
| ares-solid | Consistent body-force loads | delivered (`verification::mms`) |
| ares-solid | Typed Dirichlet conditions | delivered (`oracles::dirichlet`) |
| ares-solid | Consistent traction loads on validated nodes (ADR 0002) | delivered (`oracles::neumann`); its undocumented `D = 0` panic is the verification item `ARES-FACET-DIM-ZERO` |
| ares-solid | Moirai parallel assembly, Mnemosyne arenas, Themis placement | non-goal (veto point 3) |
| ares-operator | Athena `LinearOperator` over Leto, conjugate-gradient solve | delivered (`ares-operator` `operator` tests) |
| ares-operator | Manufactured solution at `O(h^2)`, energy identity, cantilever, Lamé cylinder | delivered at `f64` (`verification::*`); `f32` for MMS accuracy and energy (`scalar.rs`), whose rule exempts rate studies only (`scalar.rs:18-27`) |
| ares-operator | Lamé cylinder at `f32` | in scope, `M-vessel-wall`: `vessel_wall` runs at both scalars |
| ares-operator | Cantilever at `f32`, single-mesh accuracy (condition 5) | in scope, `M-typed-exchange`: test `scalar::the_cantilever_matches_beam_theory_at_f32` |
| ares-operator | Solve on Hephaestus `HostDevice` through `CsrOperator` | in scope, `M-vessel-wall`: test `the_host_device_solve_matches_the_host_solve`, and `vessel_wall` |
| ares-operator | Solve on an accelerator (wgpu, CUDA, Metal, ROCm) | non-goal (veto point 2) |
| ares-operator | Preconditioner selection, solver policy | non-goal (Athena owns it; book chapter 8) |
| ares-coupling | One-way quasi-static traction partition, conformity rejection, work identity, zero-traction null case | delivered (`coupling` tests, `f64` only) |
| ares-coupling | Cylinder under coupled pressure, interface from a Gaia boundary, work identity, rigid-translation invariance, `f32` and `f64` (ADR 0059 oracles) | in scope, `M-vessel-wall`: `vessel_wall` self-checks each, at both scalars |
| ares-coupling | `Stress`-typed traction envelope into the partition | in scope, `M-vessel-wall`: test `a_stress_envelope_loads_a_strip_through_the_partition` |
| ares-coupling | `Length`-typed displacement envelope, marshalling round trip | in scope, `M-typed-exchange`: test `a_traction_envelope_round_trips_through_export_and_advance` |
| ares-coupling | Interface stress diagnostic envelope | in scope, `M-typed-exchange`: example `typed_exchange` checks it against a closed-form stress computed without the solve |
| ares-coupling | Two-way coupling, ALE, non-conforming transfer, fluid traction extraction | non-goal (ADR 0059) |
| workspace | Book, eleven chapters; crate READMEs | book delivered; `mdbook test` and the README build are 1.0 standing conditions (`ARES-CRATE-README-DOCS`) |
| workspace | Publication of `ares-solid`, then `ares-operator` | 1.0 standing condition (`cargo publish --dry-run`, ADR 0001:82-83); the release itself needs release authority |
| workspace | Publication of `ares-coupling` | non-goal for 1.0: it depends on `harmonia` (`publish = false`), whose graduation the ADR 0001 revision of 2026-09-04 places outside this repository (`gap_audit.md`) |
| workspace | Dynamics, finite deformation, plasticity, viscoelasticity, contact, fracture, fatigue, anisotropy, composites, thermal expansion, buckling | non-goal (later phases, each with its own charter) |
| workspace | CFDrs first-consumer integration; Kwavers migration | outside this repository (atlas epic; Phase 1) |

## Hypotheses

- **Host device.** `CsrOperator` needs `SparseOperatorOps<D, T>` and
  `DenseVectorOps<D, T>` (athena 041afcb `operator.rs:67-72`); hephaestus-host
  implements both for `HostDevice` (hephaestus 742e5fb9 `sparse.rs:71`,
  `dense_vector.rs:124`), both on `hephaestus-core` 0.28.0. Not compiled: the
  `f32` and `f64` scalar bounds are unchecked. `HostDevice` needs no hardware,
  so the check never skips. Verification step: `ARES-SPIKE-HOSTDEVICE-CSR`
  compiles this composition before adding public API.
- **Gaia input.** `TetrahedralMeshBuilder` builds a positively oriented
  `IndexedMesh`. Not yet shown: whether it lends `[T; 3]` coordinates and
  `[usize; 4]` connectivity or needs a copy into `SimplexMesh` slices.
  Verification step: `ARES-SPIKE-GAIA-ACCESSOR` reads its accessors.
- **Typed envelope.** `FieldEnvelope` holds one `Quantity` per cell of a
  Cartesian `GridGeometry` (harmonia 01ebc48 `envelope.rs:42`), so a curved
  Gaia interface cannot carry one and `vessel_wall` keeps the flat partition
  slice. Not yet shown: how a `D`-component traction on a straight strip maps
  onto grid cells under ADR 0059's node-major order (0059:116-117).
  `Substep::new` is `pub(crate)`, so the round trip may need an upstream
  constructor. Verification steps: `ARES-SPIKE-STRIP-ENVELOPE` fixes the
  layout before any API; `ARES-SPIKE-SUBSTEP-ADVANCE` decides the constructor.

## Alternatives rejected

- **1.0 at the current surface.** The crate README would still claim principal
  stresses, and nothing would show a path through Gaia or a device.
- **Device solve as its own later milestone.** It leaves the largest risk
  second and the first milestone short of a named seam.
- **Device-resident element kernels instead of CSR.** That puts per-backend
  kernels in a balance package, with no oracle until a host CSR path exists.

## Consequences

- Open items carry `milestone:`; each feature item names its exit check.
- A milestone closes when its exit checks pass at a fetched `main` commit.
  Releases stay behind release authority.
- Removing, weakening or moving an exit check, or splitting a milestone, needs
  a revision note here and a judge verdict.

## Veto points

Scope that no source determines, recorded with the recommended option:

1. The Hephaestus route is assembled CSR through `CsrOperator`, not
   device-resident element kernels. ADR 0057 requires the seam, not a route.
   Book chapter 6 declines an assembled matrix and chapter 8 defers one; this
   builds it for the device path only, and the host path stays matrix-free.
2. 1.0 crosses the Hephaestus seam on `HostDevice` only. Accelerator runs
   test per-device kernels, which `hephaestus-conformance` owns.
3. The optional Moirai, Mnemosyne and Themis substrate (ADR 0057:71-72) is not
   1.0 scope.

## Overturning evidence

- `CsrOperator` does not compose with `HostDevice` for an Ares operator. The
  route is revised here before the first milestone closes.
- No `FieldEnvelope` layout carries a strip's traction. The first milestone's
  typed-envelope check is then revised here.
- Atlas ADR 0059 is rejected or reshaped.
- `harmonia` becomes publishable; publishing `ares-coupling` re-enters scope.
