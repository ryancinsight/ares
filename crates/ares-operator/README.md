# ares-operator

The [Athena](https://github.com/ryancinsight/athena) linear-operator seam for
the [Ares](https://github.com/ryancinsight/ares) solid momentum balance.

Ares assembles `f = K u` matrix-free over a mesh. Athena's Krylov solvers take
a `LinearOperator` that answers exactly that question. This crate is the join —
it owns no physics, no discretisation, and no solver policy.

## Why it is a separate crate

`ares` is `#![no_std]` and depends on nothing but vocabulary crates. Athena's
operator trait fixes the error type to its backend's, so an implementation must
name a concrete backend, and the only host backend links `std` through `leto`.
Implementing the seam inside `ares` would push that dependency into the domain
core; gating it behind a cargo feature would make the shipped configuration the
one CI does not build by default. See
[ADR 0001](https://github.com/ryancinsight/ares/blob/main/docs/adr/0001-athena-seam-as-a-separate-crate.md).

## Usage

```toml
[dependencies]
ares-operator = "0.3.0"
```

The example solves a unit square of two triangles, held along `x = 0` and
pulled along `x = 1`. It also uses `ares` (`ares-solid`), `athena-core`,
`athena-leto`, `leto`, `aequitas` and `proteus-mat`.

```rust
use aequitas::systems::si::quantities::{Dimensionless, Pressure};
use ares::{DirichletConditions, PrescribedDisplacement, SimplexMesh};
use ares_operator::ConstrainedStiffness;
use athena_core::{Cg, CgWorkspace, ConvergencePolicy, Identity};
use athena_leto::LetoBackend;
use leto::Array1;
use proteus::IsotropicModuli;

let nodes = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
let cells = [[0, 1, 2], [0, 2, 3]];
// Nodes 0 and 3 lie on `x = 0`; both of their components are held.
let prescribed = [
    PrescribedDisplacement::new(0, 0, 0.0),
    PrescribedDisplacement::new(0, 1, 0.0),
    PrescribedDisplacement::new(3, 0, 0.0),
    PrescribedDisplacement::new(3, 1, 0.0),
];
let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid mesh");
let conditions = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");
let moduli = IsotropicModuli::from_young_poisson(
    Pressure::from_base(200e9),
    Dimensionless::from_base(0.3),
)
.expect("inside the positive-definite domain");
let operator = ConstrainedStiffness::new(mesh, moduli, conditions);

// 1 MPa in +x on the unit edge between nodes 1 and 2, split evenly.
let mut external = [0.0; 8];
external[2] = 0.5e6;
external[4] = 0.5e6;
let mut load = [0.0; 8];
operator.load(&external, &mut load).expect("one entry per degree of freedom");

let backend = LetoBackend::<f64>::default();
let right_hand_side = Array1::from_shape_vec([8], load.to_vec()).expect("valid vector");
let mut solution = Array1::zeros([8]);
let mut workspace = CgWorkspace::new(&backend, 8).expect("workspace");
let policy = ConvergencePolicy::<f64>::new(1e-18, 1e-12, 200).expect("valid policy");
let report = Cg::<LetoBackend<f64>>::solve_into(
    &backend, &operator, &Identity, &right_hand_side,
    &mut solution, &mut workspace, policy,
)
.expect("well-shaped system");
assert!(report.converged());

// The held components stay exactly zero and the pulled edge moves in +x.
let displacement = solution.as_slice().expect("contiguous");
assert_eq!(displacement[..2], [0.0, 0.0]);
assert_eq!(displacement[6..], [0.0, 0.0]);
assert!(displacement[2] > 0.0 && displacement[4] > 0.0);
```

## Why the operator carries no fallible path of its own

`LinearOperator::apply` must return the backend's error type. For the host
backend that is `LetoBackendError`, a `#[non_exhaustive]` enum with no variant
for a degenerate element (`athena-leto` at `2be0542`, `src/error.rs:4-50`).
`SimplexMesh::try_new` establishes that every cell integrates, and the
Dirichlet conditions are validated against the mesh, so by the time an
operator exists the only failures left are shape mismatches, which
`LetoBackendError::LengthMismatch` names exactly.

## Residual floor

`ConstrainedStiffness::residual_floor` bounds the rounding error of evaluating
`b - A u` at the displacement that solves `A u = b`, for use as an absolute
threshold in Athena's convergence policy. It evaluates the element action's
own chain of operations on absolute values, one pass over the cells, and its
operation count mirrors `ares::stiffness_action`.

## Verification

`tests/operator.rs` checks **relay fidelity** — what Athena receives through
`apply` is bitwise identical to what `SimplexMesh::constrained_action`
produces, so the adapter cannot quietly add or lose anything — and
**solvability**: conjugate gradients converges on the operator and the field it
returns satisfies the system independently of the solver's own residual. It
also checks the prescribed-displacement path and hand-evaluated residual
floors.

`tests/verification` holds the end-to-end oracles of the solid balance, which
need a solve and so live here rather than in `ares`: a manufactured solution
with `O(h^2)` refinement, cantilever tip deflection, strain energy against
external work, and the Lame thick-walled cylinder. The manufactured solution
and the energy identity also run at `f32`; the others run at `f64`.

## Licence

MIT OR Apache-2.0.
