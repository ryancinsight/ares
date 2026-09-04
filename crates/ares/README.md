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
ares = { package = "ares-solid", version = "0.1.0" }
```

## Scope

Phase 0 is **small-strain linear elastostatics on an unstructured mesh**
(atlas ADR 0057):

- kinematics — displacement gradient and small-strain tensor;
- stress — Cauchy stress, invariants, von Mises, principal stresses;
- constitutive coupling — isotropic Hooke over `proteus::IsotropicModuli`;
- balance — static equilibrium residual;
- discretisation — continuous-Galerkin linear simplices on Gaia meshes;
- boundary conditions — Dirichlet and Neumann as typed conditions;
- assembly and solve — through Athena, backend-neutral.

Not in Phase 0: plasticity, viscoelasticity, hyperelasticity, finite
deformation, contact, dynamics, fracture, fatigue, anisotropy, buckling. Each
is a later phase with its own charter, and none is scaffolded — a module for a
capability that does not exist is a placeholder.

## Verification

Phase 0 is verified against analytical oracles rather than a reference
implementation, because none exists to difference against: the FEM patch test
and rigid-body motion exactly, then the Lamé thick-walled cylinder, cantilever
tip deflection, a manufactured solution, `O(h^2)` convergence, and strain
energy against external work. Every oracle runs at `f32` and `f64`.

## Substrate

`aequitas` quantities, `eunomia` scalars, `leto` arrays, `proteus` closure,
`gaia` mesh, `athena` solve. No `nalgebra`, `ndarray`, `rayon`, or
`num-traits`: each duplicates a capability the stack owns first-party, and
`deny.toml` enforces that rather than review.

## License

MIT OR Apache-2.0.
