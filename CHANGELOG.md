# Changelog

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
