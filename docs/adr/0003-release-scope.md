# 0003 — Release scope for 1.0

- **Status:** Proposed
- **Date:** 2026-11-14
- **Driver:** board item `ARES-RELEASE-LADDER`, delivered by PR #85
- **Class:** `[patch]` — records scope and the milestone ladder; changes no
  code.

## Context

Ares has no release tag; the workspace manifest is at 0.3.0. Its scope comes
from atlas ADR 0057 (Phase 0 charter, Accepted), atlas ADR 0059 (fluid-structure
coupling Phase 0, Proposed), atlas ADR 0055 (Proteus closes, Ares balances),
ADRs 0001 and 0002 here, the README files, the book, and the board. No source
states what 1.0 is, so this record does.

## Decision

1.0 is the complete Phase 0 vertical slice of ADR 0057 plus the one-way
coupling of ADR 0059. A user meshes a solid in Gaia, closes it with Proteus
moduli, loads it directly or through a Harmonia partition, solves through
Athena on a host or a Hephaestus device, and reads displacement and stress
measures. Typed field envelopes carry the exchange. Later phases are not
1.0 scope.

The ladder in `backlog.md` has three milestones: `M-vessel-wall` (0.3.0),
`M-device-solve` (0.4.0) and `M-typed-exchange` (1.0.0). With no release tag,
the first block takes the manifest version. The first milestone runs one real
path through every seam the scope names: Gaia mesh, Proteus closure, Athena
solve, Harmonia partition. Device solve comes next because its feasibility is
the largest open risk. Typed exchange comes last, since it rewrites the
coupling API that 1.0 freezes and ADR 0059 is still Proposed.

## Capability areas

| Crate | Area | State |
| --- | --- | --- |
| ares-solid | Displacement gradient, small strain, symmetric tensor algebra | delivered (`oracles::kinematics`) |
| ares-solid | Cauchy stress, mean and deviatoric stress, typed `J2`, von Mises | delivered (`oracles::constitutive`) |
| ares-solid | Principal stresses | in scope, `M-vessel-wall`: test `constitutive::principal_stresses_of_pure_shear_are_plus_and_minus_the_shear_stress`; the crate README claims it at base |
| ares-solid | Isotropic Hooke over `proteus` moduli; no material constants | delivered (`constitutive::the_crate_names_no_material`) |
| ares-solid | Equilibrium residual, matrix-free internal forces | delivered (`oracles::patch_test`, `oracles::assembled_operator`) |
| ares-solid | Linear simplices in 2-D and 3-D, exact one-point integration | delivered (`oracles::element_*`); for linear simplices, ADR 0057's isoparametric mapping and Gauss quadrature reduce to this |
| ares-solid | Higher-order elements, general quadrature rules | non-goal |
| ares-solid | Mesh validation over a borrowed `SimplexMesh` | delivered (`oracles::mesh_validation`) |
| ares-solid | Gaia tetrahedral mesh as solve input | in scope, `M-vessel-wall`: example `vessel_wall` |
| ares-solid | Consistent body-force loads | delivered (`verification::mms`) |
| ares-solid | Typed Dirichlet conditions | delivered (`oracles::dirichlet`) |
| ares-solid | Consistent traction loads on validated nodes (ADR 0002) | delivered (`oracles::neumann`); its undocumented `D = 0` panic is the verification item `ARES-FACET-DIM-ZERO` (`compile_fail` doctests) |
| ares-solid | Moirai parallel assembly, Mnemosyne arenas, Themis placement | non-goal (ADR 0057 lists them as optional; no measured need) |
| ares-operator | Athena `LinearOperator` over Leto, conjugate-gradient solve | delivered (`ares-operator` `operator` tests) |
| ares-operator | Lamé cylinder, cantilever, manufactured solution at `O(h^2)`, energy identity | delivered (`verification::*`); `f32` for accuracy and identity, `f64` for rate studies |
| ares-operator | Solve on a Hephaestus device | in scope, `M-device-solve`: example `device_cylinder`, test `the_device_solve_matches_the_host_solve` |
| ares-operator | Preconditioner selection, solver policy | non-goal (Athena owns it; book chapter 8) |
| ares-coupling | One-way quasi-static traction partition, conformity rejection, work identity, zero-traction null case | delivered (`coupling` tests, `f64` only) |
| ares-coupling | Cylinder under coupled pressure, interface from a Gaia boundary, rigid-translation invariance, `f32` and `f64` (ADR 0059 oracles) | in scope, `M-vessel-wall`: example `vessel_wall` self-checks each, at both scalars |
| ares-coupling | `Stress`-typed traction, `Length`-typed displacement, marshalling round trip | in scope, `M-typed-exchange`: test `a_traction_envelope_round_trips_through_export_and_advance` |
| ares-coupling | Interface stress diagnostic envelope | in scope, `M-typed-exchange`: example `typed_exchange` checks it against Hooke of the solved strain |
| ares-coupling | Two-way coupling, ALE, non-conforming transfer, fluid traction extraction | non-goal (ADR 0059) |
| workspace | Book, eleven chapters | delivered; `mdbook test` is a 1.0 standing condition |
| workspace | Publication of `ares-solid` | 1.0 standing condition (`cargo publish --dry-run`); the release itself needs release authority (atlas promotion step A9) |
| workspace | Publication of `ares-coupling` | non-goal for 1.0: it depends on `harmonia`, which is `publish = false` (`gap_audit.md`) |
| workspace | Publication of `ares-operator` | non-goal for 1.0 (veto point 5); whether its `athena` and `leto` providers are on crates.io was not looked up |
| workspace | Dynamics, finite deformation, plasticity, viscoelasticity, contact, fracture, fatigue, anisotropy, composites, thermal expansion, buckling | non-goal (later phases, each with its own charter) |
| workspace | CFDrs first-consumer integration; Kwavers migration | outside this repository (atlas epic; Phase 1) |

## Hypotheses

- **Device solve.** `athena-hephaestus` offers `CsrOperator`, a device-resident
  sparse operator, and no matrix-free path. Ares assembles matrix-free on the
  host. The recommended route assembles the constrained stiffness into CSR, so
  the host matrix-free solve is the oracle for the device solve. Not yet shown:
  that Athena accepts that operator for this problem, and that the gate host
  has a Hephaestus device. Verification step: the first item of
  `M-device-solve` probes both before any public API is added.
- **Gaia input.** Gaia's `TetrahedralMeshBuilder` builds an `IndexedMesh` with
  positive cell orientation. Not yet shown: whether it lends contiguous
  `[T; 3]` coordinates and `[usize; 4]` connectivity, or whether a copy into
  `SimplexMesh` slices is needed. Verification step: read its accessors at the
  first `M-vessel-wall` item.
- **Typed exchange.** Harmonia's `FieldEnvelope::try_new` is public, while
  `Substep::new` is `pub(crate)`. The atlas promotion item records the
  inherent-method workaround this forces. A round-trip test through `advance`
  may need an upstream constructor.

## Alternatives rejected

- **1.0 at the current surface.** The crate README would still claim principal
  stresses, and nothing would show a path through Gaia or a device.
- **Typed exchange first.** It rests on a Proposed ADR, and it would leave the
  device-solve risk until last.
- **Device-resident element kernels instead of CSR.** That puts per-backend
  kernels in a balance package, and there is no oracle until a host CSR path
  exists.

## Consequences

- Open items carry `milestone:`; `ARES-README-CRATE-TABLE` does not, because
  PR #84 deletes its block. Each feature item names the exit check it serves.
- A milestone closes when its exit checks pass at a fetched `main` commit.
  Releases stay behind release authority.
- Removing, weakening or moving an exit check, or splitting a milestone, needs
  a revision note here and a judge verdict.

## Veto points

Scope that no source determines, recorded with the recommended option:

1. Gaia mesh input belongs in the first milestone (0.3.0), not a later one.
2. The Hephaestus device solve is 1.0 scope, through assembled CSR. ADR 0057
   requires backend neutrality; book chapter 8 defers assembled-matrix work.
3. Typed envelope exchange is 1.0 scope while atlas ADR 0059 is Proposed.
4. The optional Moirai, Mnemosyne and Themis substrate is not 1.0 scope.
5. 1.0 publishes `ares-solid` only. `ares-operator` and `ares-coupling` then
   need `publish = false` in their manifests, which an item files.

## Overturning evidence

- The device probe shows that `athena-hephaestus` cannot serve an assembled
  Ares operator. The device milestone then moves past 1.0, and this record is
  revised.
- Atlas ADR 0059 is rejected or reshaped.
- `harmonia` becomes publishable. Publishing `ares-coupling` then re-enters
  scope.
