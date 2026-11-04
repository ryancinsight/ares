# Changelog

## [0.3.0] - Unreleased

### Breaking

- [major] `TractionBoundary` retains its validated coordinates by immutable
  borrow. Replace `boundary.add_consistent_loads(&nodes, &mut loads)` with
  `boundary.add_consistent_loads(&mut loads)`. Keep the coordinates alive and
  immutable while using the boundary; after changing geometry, construct a
  new boundary with `try_new`. The workspace components move together to
  0.3.0; external provider requirements and lock records stay unchanged.
- [major] `facet_measure` and `TractionBoundary::try_new` reject `D = 0` at
  compilation, including empty boundaries. Positive dimensions remain generic.

### Fixed

- A validated boundary cannot assemble onto a foreign coordinate array. The
  former replacement-array argument allowed an out-of-bounds panic with fewer
  nodes and silently zeroed traction on collapsed replacement geometry. Load
  shape errors now refer to the retained domain and leave all loads unchanged.

### Release status

- Source validation does not establish registry compatibility. The normalized
  0.2.0 baseline fails its generic material oracle to compile: published
  Proteus 0.1.1 uses Eunomia 0.10 while Ares uses 0.12. Compatible provider
  publication/adoption and complete component package qualification remain
  prerequisites; see the [active board](backlog.md).

## [0.2.0] - Unreleased

### Changed

- [minor] `CauchyStress::second_deviatoric_invariant` returns the new
  `PressureSquared` quantity instead of a raw scalar: `J2` is `Pressure ×
  Pressure` through the aequitas dimension algebra, and `von_mises` takes
  `sqrt(3 J2)` through the type-level root, so the dimension relation is
  checked at compile time. Closed-form oracles pin `sigma^2/3` (uniaxial
  stress) and `tau^2` (pure shear).

## [Unreleased]

### Added

- Phase 0 of the solid momentum balance (atlas ADR 0057): kinematics, isotropic
  Hooke over a Proteus closure, linear simplex elements, matrix-free global
  assembly, typed Dirichlet and Neumann conditions, and the Athena
  linear-operator seam.

  The patch test — a constant-strain field reproduced on an arbitrary distorted
  patch to machine precision — passes in 2-D and 3-D, under pure shear and pure
  dilation, at `f32` and `f64`. Element stiffness columns additionally agree
  with hand computation through the Voigt route, which shares no code with the
  tensor formulation the implementation uses.

- Consistent nodal loads for a body force, and end-to-end verification against
  analytical oracles: a manufactured solution recovering second-order
  convergence, Lame's thick-walled cylinder, cantilever tip deflection against
  beam theory, and the strain-energy-equals-external-work identity. The
  accuracy and identity oracles run at `f32` as well as `f64`; the
  convergence-rate studies stay at `f64`, because `f32` reaches its precision
  floor before the study leaves the asymptotic regime.

- One-way fluid-to-solid coupling (atlas ADR 0059) in `ares-coupling`: a
  `StructuralInterface` carrying facet-major traction and node-major
  displacement, and a `StructuralPartition` implementing Harmonia's
  `Partition`. Interface work is conserved exactly — the work the fluid does
  equals twice the stored strain energy — which holds because the nodal load is
  the consistent one; a lumped load with the same resultant breaks it, and the
  suite was mutation-checked against exactly that substitution.

  Non-conforming interfaces are rejected with a typed error rather than
  transferred approximately.

### Changed

- The repository is a workspace: `crates/ares` (`ares-solid`) is the `no_std`,
  allocation-free domain core, and `crates/ares-operator` is the Athena operator
  seam, which links `std` through `leto`. See
  [ADR 0001](docs/adr/0001-athena-seam-as-a-separate-crate.md).

- Repository scaffold at the stack lint and gate floor: pinned toolchain,
  committed nextest budgets, `deny.toml` carrying the atlas ADR 0055 substrate
  prohibitions, pedantic clippy with `unwrap_used`, `indexing_slicing`, and
  `arithmetic_side_effects` denied, and `#![forbid(unsafe_code)]` with
  `#![deny(missing_docs)]`.

  No physics yet. The scaffold passes the full gate before any is added, so a
  later failure is attributable to the physics rather than to the floor.
