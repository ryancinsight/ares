# 8. The solve, through Athena

The assembled operator and the load vector give a linear system. Ares does not
solve it: [Athena](https://github.com/ryancinsight/athena) owns solver policy,
and `ares-athena` is the seam between them.

## Why conjugate gradients

The constrained operator is symmetric and positive definite (chapter 7), which
is exactly the class conjugate gradients is for. CG needs only the *action*
`K u`, never the matrix — which is what allowed chapter 6 to skip assembling
one.

## The seam

```rust,ignore
use ares::{DirichletConditions, SimplexMesh};
use ares_athena::ConstrainedStiffness;
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

`ConstrainedStiffness` implements Athena's `LinearOperator`. It borrows the
mesh, conditions, and material, and owns only the scratch buffer the
constrained action needs — allocated once, because a Krylov solve applies the
operator every iteration.

## Two facts about Athena's trait that shaped everything upstream

Athena's seam is:

```text
trait LinearOperator<B: KrylovBackend> {
    fn dimension(&self) -> usize;
    fn apply(&self, backend: &B, input: B::View<'_>, output: B::ViewMut<'_>)
        -> Result<(), B::Error>;
}
```

**The views are backend-associated types.** An implementation generic over `B`
receives an opaque `B::View<'_>` with no method for reading element data. So a
generic implementation cannot exist — the seam is implementable only against a
*named* backend. The only host backend links `std`, which is why the seam is a
separate crate and the domain core stays `no_std`.

**The error is fixed to `B::Error`.** An implementation cannot introduce a
failure mode the backend does not already name, and `LetoBackendError` has no
variant for "this element is degenerate".

That second constraint is why chapter 6 validates as thoroughly as it does. It
was not worked around: `SimplexMesh::try_new` establishes that every cell
integrates, and the conditions are validated against the mesh, so by the time
an operator exists the only reachable failure is a shape mismatch — which
`LengthMismatch` names exactly. A downstream constraint decided an upstream
design.

## Reading the report

`Ok(_)` does **not** mean converged. Athena returns numerical termination —
breakdown, non-positive curvature, budget exhaustion — value-semantically in
the report, and `Err` only for dimension mismatches and backend failures.

`SolveReport` is `#[must_use]` precisely to force the check. An unconverged
displacement field is smooth and plausible; nothing about it looks wrong.

## Preconditioning

The examples use `Identity` — no preconditioning. That is fine for the small
problems the verification suite runs and will not be for large ones: the
condition number of an elasticity operator grows as the mesh refines, and CG's
iteration count grows with its square root.

Athena provides Jacobi, incomplete LU, and SOR preconditioners, all of which
want the assembled matrix that chapter 6 declines to build. Reconciling those
two facts — matrix-free assembly against matrix-requiring preconditioners — is
real work and is not part of Phase 0. The honest statement is that Phase 0
solves problems small enough not to need it, and that a large problem will need
this question answered rather than ignored.
