//! Executable evidence for the Athena operator seam (atlas ADR 0057, A5).
//!
//! Two things are being checked, and they are different. First, that the seam
//! is *faithful*: what Athena receives through `LinearOperator::apply` is
//! exactly what `SimplexMesh::constrained_action` produces, bit for bit. A
//! transposed view, a stride, or an off-by-one in the adapter would otherwise
//! surface only as a solver that converges to the wrong field.
//!
//! Second, that the operator is *solvable*: Athena's conjugate gradients
//! converges on it, and the field it returns satisfies the system it was given.
//! The analytical oracles — Lame, cantilever, manufactured solutions,
//! convergence order — belong to A6; this asks only whether the seam carries a
//! solve at all.

#![expect(
    clippy::float_cmp,
    reason = "both exact comparisons are exact by derivation rather than by luck. The relay check compares two routes through the same arithmetic, so any difference is an adapter defect and not rounding. The constrained displacements are exactly zero by induction over the iteration: the initial guess is zero, so the constrained residual starts at b_c - u_c = 0; the operator's constrained rows are the identity, so (A p)_c = p_c; and every update is therefore a multiply-add on exact zeros, which stays exactly zero for as many iterations as the solver runs."
)]

use aequitas::systems::si::quantities::{Dimensionless, Pressure};
use ares::{
    DirichletConditions, PrescribedDisplacement, SimplexMesh, TractionBoundary, TractionFacet,
};
use ares_operator::ConstrainedStiffness;
use athena_core::{Cg, CgWorkspace, ConvergencePolicy, Identity, KrylovBackend, LinearOperator};
use athena_leto::LetoBackend;
use eunomia::RealField;
use leto::Array1;
use proteus::IsotropicModuli;

fn moduli<T: RealField>(young: f64, poisson: f64) -> IsotropicModuli<T> {
    IsotropicModuli::from_young_poisson(
        Pressure::from_base(T::from_f64(young)),
        Dimensionless::from_base(T::from_f64(poisson)),
    )
    .expect("inside the positive-definite domain")
}

/// Nodes, cells, and conditions for the shared fixture below.
type Fixture = (
    [[f64; 2]; 4],
    [[usize; 3]; 2],
    [PrescribedDisplacement<f64>; 4],
);

/// A unit square split into two triangles, fixed along `x = 0`.
fn fixture() -> Fixture {
    let nodes = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let cells = [[0, 1, 2], [0, 2, 3]];
    // Nodes 0 and 3 are the x = 0 edge; both components held.
    let conditions = [
        PrescribedDisplacement::new(0, 0, 0.0),
        PrescribedDisplacement::new(0, 1, 0.0),
        PrescribedDisplacement::new(3, 0, 0.0),
        PrescribedDisplacement::new(3, 1, 0.0),
    ];
    (nodes, cells, conditions)
}

#[test]
fn the_operator_reports_the_mesh_dimension() {
    let (nodes, cells, prescribed) = fixture();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid mesh");
    let bc = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");
    let operator = ConstrainedStiffness::new(mesh, moduli::<f64>(200e9, 0.3), bc);
    assert_eq!(operator.dimension(), 8);
    assert_eq!(operator.dimension(), mesh.degrees_of_freedom());
}

#[test]
fn the_seam_relays_the_constrained_action_exactly() {
    // Faithfulness. The adapter must add nothing and lose nothing, so the two
    // routes are compared for bitwise equality rather than within a tolerance:
    // any difference at all is an adapter defect, not rounding.
    let (nodes, cells, prescribed) = fixture();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid mesh");
    let bc = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");
    let material = moduli::<f64>(200e9, 0.3);
    let operator = ConstrainedStiffness::new(mesh, material, bc);

    let values: Vec<f64> = (0..8)
        .map(|i| 1e-4 * (f64::from(i * 5 % 7) - 3.0))
        .collect();
    let backend = LetoBackend::<f64>::default();
    let input = Array1::from_shape_vec([8], values.clone()).expect("valid vector");
    let mut through_athena = Array1::zeros([8]);
    operator
        .apply(
            &backend,
            backend.view(&input),
            backend.view_mut(&mut through_athena),
        )
        .expect("well-shaped");

    let mut direct = [0.0_f64; 8];
    let mut scratch = [0.0_f64; 8];
    mesh.constrained_action(&material, &bc, &values, &mut direct, &mut scratch)
        .expect("well-shaped");

    assert_eq!(
        through_athena.as_slice().expect("contiguous"),
        direct.as_slice(),
        "the adapter is not relaying the constrained action unchanged"
    );
}

#[test]
fn conjugate_gradients_converges_on_the_constrained_operator() {
    // Solvability. A traction on the free edge, the fixed edge held, and the
    // system solved through Athena's PCG with no preconditioner.
    let (nodes, cells, prescribed) = fixture();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid mesh");
    let bc = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");
    let material = moduli::<f64>(200e9, 0.3);
    let operator = ConstrainedStiffness::new(mesh, material, bc);

    // Pull the x = 1 edge in +x.
    let facets = [TractionFacet::new([1, 2], [1.0e6_f64, 0.0])];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("valid facet");
    let mut external = [0.0_f64; 8];
    boundary
        .add_consistent_loads(&nodes, &mut external)
        .expect("well-shaped");

    let mut load = [0.0_f64; 8];
    operator.load(&external, &mut load).expect("well-shaped");

    let backend = LetoBackend::<f64>::default();
    let right_hand_side = Array1::from_shape_vec([8], load.to_vec()).expect("valid vector");
    let mut solution = Array1::zeros([8]);
    let mut workspace = CgWorkspace::new(&backend, 8).expect("workspace");
    let policy = ConvergencePolicy::<f64>::new(1e-18, 1e-12, 200).expect("valid policy");

    let report = Cg::<LetoBackend<f64>>::solve_into(
        &backend,
        &operator,
        &Identity,
        &right_hand_side,
        &mut solution,
        &mut workspace,
        policy,
    )
    .expect("well-shaped system");
    assert!(
        report.converged(),
        "conjugate gradients did not converge: {:?} after {} iterations, residual {:.3e}",
        report.termination,
        report.iterations,
        report.final_residual_norm
    );

    let displacement = solution.as_slice().expect("contiguous");

    // The solution satisfies the system it was handed. Checked independently
    // of the solver's own residual, which is the quantity it minimised and so
    // cannot be evidence about the operator.
    let mut image = [0.0_f64; 8];
    let mut scratch = [0.0_f64; 8];
    mesh.constrained_action(&material, &bc, displacement, &mut image, &mut scratch)
        .expect("well-shaped");
    let load_scale = load.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    for (dof, (got, want)) in image.iter().zip(load.iter()).enumerate() {
        assert!(
            (got - want).abs() <= load_scale * 1e-9,
            "dof {dof}: K u = {got} but b = {want}"
        );
    }

    // The constrained nodes did not move, and the free ones did. Without the
    // second half every assertion above is satisfied by the zero field.
    for condition in &prescribed {
        assert_eq!(displacement[condition.degree_of_freedom::<2>()], 0.0);
    }
    assert!(
        displacement[2] > 0.0 && displacement[4] > 0.0,
        "the loaded edge did not move in the direction it was pulled: {displacement:?}"
    );
}

#[test]
fn a_prescribed_displacement_is_reproduced_by_the_solve() {
    // The non-zero Dirichlet path end to end: prescribe a displacement, solve
    // with no external load, and the solve must return the prescribed values
    // at the constrained nodes. This is what the `K u_g` term in the load
    // exists for, now checked through the solver rather than at the load.
    let nodes = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    let cells = [[0, 1, 2], [0, 2, 3]];
    let stretch = 3.0e-4_f64;
    let prescribed = [
        PrescribedDisplacement::new(0, 0, 0.0),
        PrescribedDisplacement::new(0, 1, 0.0),
        PrescribedDisplacement::new(1, 0, stretch),
        PrescribedDisplacement::new(3, 0, 0.0),
        PrescribedDisplacement::new(3, 1, 0.0),
    ];
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid mesh");
    let bc = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");
    let operator = ConstrainedStiffness::new(mesh, moduli::<f64>(200e9, 0.3), bc);

    let mut load = [0.0_f64; 8];
    operator
        .load(&[0.0_f64; 8], &mut load)
        .expect("well-shaped");

    let backend = LetoBackend::<f64>::default();
    let right_hand_side = Array1::from_shape_vec([8], load.to_vec()).expect("valid vector");
    let mut solution = Array1::zeros([8]);
    let mut workspace = CgWorkspace::new(&backend, 8).expect("workspace");
    let policy = ConvergencePolicy::<f64>::new(1e-18, 1e-12, 200).expect("valid policy");
    let report = Cg::<LetoBackend<f64>>::solve_into(
        &backend,
        &operator,
        &Identity,
        &right_hand_side,
        &mut solution,
        &mut workspace,
        policy,
    )
    .expect("well-shaped system");
    assert!(report.converged(), "{:?}", report.termination);

    let displacement = solution.as_slice().expect("contiguous");
    for condition in &prescribed {
        let dof = condition.degree_of_freedom::<2>();
        assert!(
            (displacement[dof] - condition.value()).abs() <= stretch * 1e-9,
            "dof {dof} solved to {} but was prescribed {}",
            displacement[dof],
            condition.value()
        );
    }
    // Node 2 is free and shares a cell with the stretched node, so it must
    // have moved: a solve that merely echoed the prescription would not.
    assert!(
        displacement[4].abs() > stretch * 1e-3,
        "the prescribed motion did not propagate to the free nodes"
    );
}
