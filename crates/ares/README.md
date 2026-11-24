# Ares

Solid momentum balance for the [Atlas](https://github.com/ryancinsight/atlas)
stack: kinematics, stress measures, equilibrium, and the boundary conditions
that close them.

Ares owns no material data. Constitutive closure belongs to
[Proteus](https://github.com/ryancinsight/proteus) — Proteus closes, Ares
balances (atlas ADR 0055).

## Installation

The crates.io name `ares` belongs to an unrelated third-party crate, so this
publishes as `ares-solid`. The import path stays `ares` via `[lib] name`, so
rename the dependency and no `use ares::…` changes:

```toml
[dependencies]
ares = { package = "ares-solid", version = "0.3.0" }
```

## Example

A unit right triangle stretched by `u_x = e x`, `e = 2^-10`. Young's modulus
`5/2` and Poisson's ratio `1/4` give `lambda = mu = 1`, so the stress is
`diag(3e, e)` and each nodal force is the cell measure `1/2` times that stress
applied to the node's shape gradient. All values are dyadic, so the assertion
is exact. Two dimensions means plane strain (`docs/book/balance.md`). The
example also uses `aequitas` and `proteus-mat`, dependencies of `ares-solid`.

```rust
use aequitas::systems::si::quantities::{Dimensionless, Pressure};
use ares::SimplexMesh;
use proteus::IsotropicModuli;

let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
let mesh = SimplexMesh::try_new(&nodes, &[[0, 1, 2]]).expect("one positively wound cell");
let moduli = IsotropicModuli::from_young_poisson(
    Pressure::from_base(2.5),
    Dimensionless::from_base(0.25),
)
.expect("inside the positive-definite domain");

// Degree of freedom `node * 2 + component`; only node 1 moves, in x.
let e = 1.0 / 1024.0;
let mut forces = [0.0; 6];
mesh.internal_forces(&moduli, &[0.0, 0.0, e, 0.0, 0.0, 0.0], &mut forces)
    .expect("one entry per degree of freedom");
assert_eq!(forces, [-1.5 * e, -0.5 * e, 1.5 * e, 0.0, 0.0, 0.5 * e]);
```

## Scope

Phase 0 is **small-strain linear elastostatics on an unstructured mesh**
(atlas ADR 0057):

- kinematics — small strain from a displacement gradient, symmetric tensors;
- stress — Cauchy stress, mean stress, deviator, second invariant, von Mises;
- constitutive coupling — isotropic Hooke over `proteus::IsotropicModuli`;
- discretisation — linear simplices over a borrowed `SimplexMesh`;
- assembly — matrix-free stiffness action, cell strains, body loads;
- boundary conditions — Dirichlet and Neumann as typed conditions.

Chartered for Phase 0 and not yet implemented: principal stresses, a static
equilibrium residual, and a Gaia mesh producer. The solve is in
`ares-operator`.

Not in Phase 0: plasticity, viscoelasticity, hyperelasticity, finite
deformation, contact, dynamics, fracture, fatigue, anisotropy, buckling. Each
is a later phase with its own charter, and none is scaffolded — a module for a
capability that does not exist is a placeholder.

## Boundaries

Gaia owns mesh, geometry, and proximity queries. Athena owns solver policy.
Horae owns time when dynamics arrive. Harmonia owns every coupling to another
balance domain, so Ares carries no dependency on a fluid, acoustic, or
transport package in any phase.

## Verification

Phase 0 is verified against analytical oracles rather than a reference
implementation, because none exists to difference against. This crate's
`tests/oracles` holds the patch test (constant strain on distorted 2-D and 3-D
patches), rigid-body motion, Hooke against its closed forms, hand-computed
stiffness columns, and traction resultants. The solve-level oracles live in
`ares-operator`: a manufactured solution with `O(h^2)` convergence, strain
energy against external work, cantilever tip deflection, and the Lame
thick-walled cylinder.

Scalar coverage differs by oracle. Rigid-body motion, the manufactured
solution and the energy identity run at `f32` and `f64`. The patch test runs at
`f64` on every patch and at `f32` on the distorted 2-D patch. The Lame
cylinder, the cantilever and the convergence order run at `f64` only.

## Substrate

`aequitas` quantities, `eunomia` scalars, `proteus` closure; `leto` arrays and
`athena` solvers enter through `ares-operator`. No `nalgebra`, `ndarray`,
`rayon`, or `num-traits`: each duplicates a capability the stack owns
first-party, and `deny.toml` enforces that rather than review.

## License

MIT OR Apache-2.0.
