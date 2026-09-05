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
one CI does not build by default. See [ADR 0001](../../docs/adr/0001-athena-seam-as-a-separate-crate.md).

## Usage

```toml
[dependencies]
ares-operator = "0.1.0"
```

```rust,ignore
use ares::{DirichletConditions, SimplexMesh};
use ares_operator::ConstrainedStiffness;
use athena_core::{Cg, CgWorkspace, ConvergencePolicy, Identity};
use athena_leto::LetoBackend;

let mesh = SimplexMesh::try_new(&nodes, &cells)?;
let conditions = DirichletConditions::try_new(&prescribed, mesh.node_count())?;
let operator = ConstrainedStiffness::new(mesh, moduli, conditions);

let mut load = vec![0.0; operator.dimension()];
operator.load(&external, &mut load)?;

let report = Cg::<LetoBackend<f64>>::solve_into(
    &backend, &operator, &Identity, &right_hand_side,
    &mut solution, &mut workspace, policy,
)?;
assert!(report.converged());
```

## Verification

The seam is checked for **relay fidelity** — what Athena receives through
`apply` is bitwise identical to what `SimplexMesh::constrained_action`
produces, so the adapter cannot quietly add or lose anything — and for
**solvability**, that conjugate gradients converges on the operator and the
field it returns satisfies the system independently of the solver's own
residual.

The analytical oracles for the physics — Lame, cantilever, manufactured
solutions, convergence order — live in `ares` and in atlas ADR 0057 phase A6.

## Licence

MIT OR Apache-2.0.
