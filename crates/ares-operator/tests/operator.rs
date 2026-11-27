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
    reason = "the exact comparisons are exact by derivation rather than by luck. The relay check compares two routes through the same arithmetic, so any difference is an adapter defect and not rounding. The constrained displacements are exactly zero by induction over the iteration: the initial guess is zero, so the constrained residual starts at b_c - u_c = 0; the operator's constrained rows are the identity, so (A p)_c = p_c; and every update is therefore a multiply-add on exact zeros, which stays exactly zero for as many iterations as the solver runs. The residual-floor checks are exact too: the hand-evaluated values are short dyadic rationals, so every intermediate is exact and the expected floor applies the same `gamma` and `sqrt` as the method, and doubling is exact in binary floating point."
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
        .add_consistent_loads(&mut external)
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

/// The fixture's square as an operator, holding `held` with `moduli`.
fn held_square<'a>(
    (nodes, cells, _): &'a Fixture,
    held: &'a [PrescribedDisplacement<f64>],
    moduli: IsotropicModuli<f64>,
) -> ConstrainedStiffness<'a, f64, 2, 3> {
    let mesh = SimplexMesh::try_new(nodes, cells).expect("valid mesh");
    let bc = DirichletConditions::try_new(held, mesh.node_count()).expect("valid");
    ConstrainedStiffness::new(mesh, moduli, bc)
}

/// `gamma_k` for `k = N + 2 D + 9 + 2 = 18`: nodes 0 and 2 of the square lie
/// in both cells, so two cells meet at the busiest node.
fn gamma_18<T: RealField>() -> T {
    let one = T::from_f64(1.0);
    let unit_roundoff = T::EPSILON / T::from_f64(2.0);
    let terms = T::from_f64(18.0);
    terms * unit_roundoff / (one - terms * unit_roundoff)
}

#[test]
fn the_residual_floor_is_the_hand_evaluated_absolute_chain() {
    // E = 5/2 and nu = 1/4 give mu = lambda = 1. Every intermediate is a
    // short dyadic rational (areas 1/2, shape gradients 0 and +-1), so the
    // chain is exact and the floor is compared for equality.
    //
    // Held entries enter the stiffness as zeros, so |u| is (3, 1) at node 1
    // and (1, 2) at node 2. Cell [0,1,2] has shape gradients (-1,0), (1,-1),
    // (0,1): |G| = [[3,4],[1,3]], tr|eps| = 6, |sigma| = 2|eps| + 6 I =
    // [[12,5],[5,12]], node forces (6,2.5), (8.5,8.5), (2.5,6). Cell [0,2,3]
    // has (0,-1), (1,0), (-1,1): |G| = [[1,0],[2,0]], tr|eps| = 1, |sigma| =
    // [[3,2],[2,1]], node forces (1,.5), (1.5,1), (2.5,1.5). The free rows
    // are node 1: (8.5,8.5) and node 2: (4,7); the held rows are identity
    // rows, the signed entries (-2,1) and (2,-4). Then ||A|u||^2 = 144.5 +
    // 65 + 25 = 234.5.
    let fixture = fixture();
    let operator = held_square(&fixture, &fixture.2, moduli(2.5, 0.25));
    let displacement = [-2.0, 1.0, 3.0, -1.0, -1.0, 2.0, 2.0, -4.0];

    let floor = operator.residual_floor(&displacement).expect("well-shaped");
    assert_eq!(floor, 2.0 * gamma_18::<f64>() * 234.5_f64.sqrt());
}

#[test]
fn the_residual_floor_counts_a_free_reference_node() {
    // Only node 3 is held, so node 0, the reference of both cells, moves:
    // u = (1, 0) there and zero elsewhere, so |u_0| enters |G| and every
    // |u_b| is zero.
    // Cell [0,1,2]: |G| = [[1,2],[0,0]], tr|eps| = 1, |sigma| =
    // [[3,2],[2,1]], node forces (1.5,1), (2.5,1.5), (1,.5). Cell [0,2,3]:
    // |G| = [[2,1],[0,0]], tr|eps| = 2, |sigma| = [[6,1],[1,2]], node forces
    // (.5,1), (3,.5), (3.5,1.5). The rows are node 0: (2,2), node 1:
    // (2.5,1.5), node 2: (4,1), and node 3 held at zero, so ||A|u||^2 = 33.5.
    let fixture = fixture();
    let held = [
        PrescribedDisplacement::new(3, 0, 0.0),
        PrescribedDisplacement::new(3, 1, 0.0),
    ];
    let operator = held_square(&fixture, &held, moduli(2.5, 0.25));
    let displacement = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    let floor = operator.residual_floor(&displacement).expect("well-shaped");
    assert_eq!(floor, 2.0 * gamma_18::<f64>() * 33.5_f64.sqrt());
}

#[test]
fn the_residual_floor_adds_the_reference_and_neighbour_motion() {
    // Only node 3 is held. Nodes 0 and 1 both move by (1, 0), so in cell
    // [0,1,2] the term (|u_1| + |u_0|) is (2, 0), a sum of two nonzero
    // magnitudes that a difference or a maximum would not reproduce. Cell
    // [0,1,2]: |G| = [[2,3],[0,0]], tr|eps| = 2, |sigma| = [[6,3],[3,2]],
    // node forces (3,1.5), (4.5,2.5), (1.5,1). Cell [0,2,3]: |G| =
    // [[2,1],[0,0]], tr|eps| = 2, |sigma| = [[6,1],[1,2]], node forces
    // (.5,1), (3,.5), (3.5,1.5). The rows are node 0: (3.5,2.5), node 1:
    // (4.5,2.5), node 2: (4.5,1.5), and node 3 held at zero, so ||A|u||^2 =
    // 18.5 + 26.5 + 22.5 = 67.5.
    let fixture = fixture();
    let held = [
        PrescribedDisplacement::new(3, 0, 0.0),
        PrescribedDisplacement::new(3, 1, 0.0),
    ];
    let operator = held_square(&fixture, &held, moduli(2.5, 0.25));
    let displacement = [1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0];

    let floor = operator.residual_floor(&displacement).expect("well-shaped");
    assert_eq!(floor, 2.0 * gamma_18::<f64>() * 67.5_f64.sqrt());
}

/// The floor of the fixture square's operator, in the scalar `T`, with only
/// node 3 held and nodes 0 and 1 each moved by `(1/2, 0)`.
///
/// The reference node 0 and its neighbour node 1 both move, by a magnitude
/// other than 1, so in cell `[0,1,2]` the term `|u_1| + |u_0|` is `1` where
/// `|u_1| + |u_0|^2` gives `3/4`. A magnitude of 1 hides that square, as
/// `1 * 1 = 1`, and so does a magnitude of 0. Halving the magnitude of
/// `the_residual_floor_adds_the_reference_and_neighbour_motion` scales
/// `||A|u||^2 = 67.5` by `1/4`: cell `[0,1,2]` has `|G| = [[1,3/2],[0,0]]`,
/// `tr|eps| = 1`, `|sigma| = [[3,3/2],[3/2,1]]`, node forces `(3/2,3/4)`,
/// `(9/4,5/4)`, `(3/4,1/2)`; cell `[0,2,3]` has `|G| = [[1,1/2],[0,0]]`,
/// `tr|eps| = 1`, `|sigma| = [[3,1/2],[1/2,1]]`, node forces `(1/4,1/2)`,
/// `(3/2,1/4)`, `(7/4,3/4)`. The rows are node 0: `(7/4,5/4)`, node 1:
/// `(9/4,5/4)`, node 2: `(9/4,3/4)`, so `||A|u||^2 = 135/8`, and the squared
/// reading of the term gives `235/32`. Every intermediate is a short dyadic
/// rational, so the chain is exact at `f32` and `f64`.
fn assert_the_floor_adds_a_half_magnitude_reference_motion<T: RealField + core::fmt::Debug>() {
    let zero = T::from_f64(0.0);
    let one = T::from_f64(1.0);
    let half = T::from_f64(0.5);
    let nodes = [[zero, zero], [one, zero], [one, one], [zero, one]];
    let cells = [[0, 1, 2], [0, 2, 3]];
    let held = [
        PrescribedDisplacement::new(3, 0, zero),
        PrescribedDisplacement::new(3, 1, zero),
    ];
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid mesh");
    let conditions = DirichletConditions::try_new(&held, mesh.node_count()).expect("valid");
    let operator = ConstrainedStiffness::new(mesh, moduli::<T>(2.5, 0.25), conditions);
    let mut displacement = [zero; 8];
    displacement[0] = half;
    displacement[2] = half;

    let floor = operator.residual_floor(&displacement).expect("well-shaped");
    let two = T::from_f64(2.0);
    assert_eq!(
        floor,
        two * gamma_18::<T>() * T::from_f64(16.875).sqrt(),
        "at {}",
        core::any::type_name::<T>()
    );
}

#[test]
fn the_residual_floor_adds_a_reference_motion_of_non_unit_magnitude() {
    assert_the_floor_adds_a_half_magnitude_reference_motion::<f32>();
    assert_the_floor_adds_a_half_magnitude_reference_motion::<f64>();
}

#[test]
fn the_residual_floor_takes_the_magnitude_of_a_negative_lambda() {
    // E = 1 and nu = -1/2 give mu = 1 and lambda = -1/2, and the bound uses
    // |lambda| = 1/2. The displacement is the chain test's, so |G| and |eps|
    // are as there. Cell [0,1,2]: tr|eps| = 6, |sigma| = [[9,5],[5,9]], node
    // forces (4.5,2.5), (7,7), (2.5,4.5). Cell [0,2,3]: tr|eps| = 1,
    // |sigma| = [[2.5,2],[2,.5]], forces (1,.25), (1.25,1), (2.25,1.25). The
    // free rows are node 1: (7,7) and node 2: (3.75,5.5), with the identity
    // rows (-2,1) and (2,-4): ||A|u||^2 = 98 + 44.3125 + 25 = 167.3125.
    let fixture = fixture();
    let operator = held_square(&fixture, &fixture.2, moduli(1.0, -0.5));
    let displacement = [-2.0, 1.0, 3.0, -1.0, -1.0, 2.0, 2.0, -4.0];

    let floor = operator.residual_floor(&displacement).expect("well-shaped");
    assert_eq!(floor, 2.0 * gamma_18::<f64>() * 167.3125_f64.sqrt());
}

#[test]
fn the_residual_floor_doubles_with_the_displacement() {
    // Doubling is exact in binary floating point, and every term of the bound
    // is linear in the displacement magnitudes, so the floor doubles exactly.
    let fixture = fixture();
    let operator = held_square(&fixture, &fixture.2, moduli(2.5, 0.25));
    let displacement = [0.0, 0.0, 1.5e-4, -2.0e-5, 7.0e-5, 3.0e-5, 0.0, 0.0];
    let doubled = displacement.map(|value| value * 2.0);

    let floor = operator.residual_floor(&displacement).expect("well-shaped");
    assert!(floor > 0.0, "a moving free node must have a nonzero floor");
    assert_eq!(
        operator.residual_floor(&doubled).expect("well-shaped"),
        2.0 * floor
    );
}

#[test]
fn the_residual_floor_rejects_a_misshaped_displacement() {
    let fixture = fixture();
    let error = held_square(&fixture, &fixture.2, moduli(2.5, 0.25))
        .residual_floor(&[0.0; 7])
        .expect_err("seven entries for eight degrees of freedom");
    assert!(matches!(
        error,
        athena_leto::LetoBackendError::LengthMismatch { left: 8, right: 7 }
    ));
}
