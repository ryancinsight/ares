//! Executable evidence for boundary conditions (atlas ADR 0057, phase A5).
//!
//! # What each oracle is for
//!
//! **Neumann.** A consistent nodal load is not "the traction split between the
//! nodes" — it is `integral(N_a t) dS`, and the two coincide only for linear
//! elements under a uniform traction. Two independent statics identities pin
//! it: the loads must carry the exact resultant **force** `t * A`, and the
//! exact resultant **moment** about any origin. Force alone does not
//! distinguish an equal split from an unequal one with the same total, so the
//! moment is what fixes the distribution.
//!
//! **Dirichlet.** The constrained operator must stay symmetric — conjugate
//! gradients has no convergence guarantee otherwise, and the failure is a
//! wrong answer rather than an error. The two defects this guards are both
//! silent: projecting on one side only leaves `K_cf` intact and breaks
//! symmetry, and dropping the `K u_g` term from the load silently solves the
//! problem in which every prescribed displacement is zero.
//!
//! That second one is why a non-zero prescribed value appears in these tests
//! at all. Fixed-at-zero is the common case, and it is exactly the case in
//! which the omission cannot be detected.

#![expect(
    clippy::cast_precision_loss,
    reason = "the casts turn small loop and node indices into the reals of a test field, far inside f64's exact-integer range."
)]
#![expect(
    clippy::float_cmp,
    reason = "the exact comparisons are exact by construction: a projected field's constrained entries are assigned zero and its carried entries are copied, so both are bit-identical to their sources rather than computed. The statics identities and the symmetry check, which do involve arithmetic, carry derived bounds."
)]

use aequitas::systems::si::quantities::{Dimensionless, Pressure};
use ares::{
    DirichletConditions, InvalidBoundary, InvalidConditions, PrescribedDisplacement, SimplexMesh,
    TractionBoundary, TractionFacet,
};
use eunomia::RealField;
use proteus::IsotropicModuli;

fn moduli<T: RealField>(young: f64, poisson: f64) -> IsotropicModuli<T> {
    IsotropicModuli::from_young_poisson(
        Pressure::from_base(T::from_f64(young)),
        Dimensionless::from_base(T::from_f64(poisson)),
    )
    .expect("inside the positive-definite domain")
}

fn square_patch() -> ([[f64; 2]; 5], [[usize; 3]; 4]) {
    (
        [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.37, 0.61]],
        [[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
    )
}

// ---------------------------------------------------------------------------
// Neumann: consistent nodal loads
// ---------------------------------------------------------------------------

#[test]
fn a_facet_measure_is_its_length_in_two_dimensions() {
    // A 3-4-5 edge, so the closed form is an exact integer.
    let nodes = [[1.0_f64, 2.0], [4.0, 6.0]];
    assert_eq!(ares::boundary::facet_measure(&nodes), 5.0);
}

#[test]
fn a_facet_measure_is_its_area_in_three_dimensions() {
    // A right triangle with legs 3 and 4 in a plane oblique to every axis, so
    // the Gram determinant is doing real work rather than reading off a
    // coordinate difference.
    let nodes = [[0.0_f64, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 0.0, 4.0]];
    assert_eq!(ares::boundary::facet_measure(&nodes), 6.0);
}

#[test]
fn a_traction_carries_its_exact_resultant_force() {
    // The first statics identity: sum of nodal loads == t * A.
    let nodes = [[0.0_f64, 0.0], [2.0, 0.0], [2.0, 3.0], [0.0, 3.0]];
    // The right-hand edge, length 3.
    let traction = [1.7e6_f64, -0.4e6];
    let facets = [TractionFacet::new([1, 2], traction)];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("valid facet");

    let mut loads = [0.0_f64; 8];
    boundary
        .add_consistent_loads(&nodes, &mut loads)
        .expect("well-shaped");

    for component in 0..2 {
        let total: f64 = (0..4).map(|node| loads[node * 2 + component]).sum();
        let expected = traction[component] * 3.0;
        let tolerance = expected.abs() * f64::EPSILON * 8.0;
        assert!(
            (total - expected).abs() <= tolerance,
            "component {component}: resultant {total} != {expected}"
        );
    }
}

#[test]
fn a_traction_carries_its_exact_resultant_moment() {
    // The second, independent identity. A uniform traction acts through the
    // facet centroid, so `sum_a x_a (x) f_a` must equal `A * x_centroid (x) t`.
    // The force check above passes for any distribution with the right total;
    // this one fixes the distribution.
    let nodes = [[0.0_f64, 0.0], [2.0, 0.0], [2.0, 3.0], [0.0, 3.0]];
    let traction = [1.7e6_f64, -0.4e6];
    let facets = [TractionFacet::new([1, 2], traction)];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("valid facet");

    let mut loads = [0.0_f64; 8];
    boundary
        .add_consistent_loads(&nodes, &mut loads)
        .expect("well-shaped");

    // Scalar moment about the origin: sum of (x f_y - y f_x).
    let moment: f64 = (0..4)
        .map(|node| nodes[node][0] * loads[node * 2 + 1] - nodes[node][1] * loads[node * 2])
        .sum();
    let centroid = [
        f64::midpoint(nodes[1][0], nodes[2][0]),
        f64::midpoint(nodes[1][1], nodes[2][1]),
    ];
    let expected = 3.0 * (centroid[0] * traction[1] - centroid[1] * traction[0]);
    let tolerance = expected.abs() * f64::EPSILON * 16.0;
    assert!(
        (moment - expected).abs() <= tolerance,
        "resultant moment {moment} != {expected}"
    );
}

#[test]
fn a_three_dimensional_traction_carries_its_exact_resultant() {
    let nodes = [
        [0.0_f64, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, 0.0, 4.0],
        [1.0, 5.0, 1.0],
    ];
    let traction = [0.0_f64, 2.5e6, 0.0];
    let facets = [TractionFacet::new([0, 1, 2], traction)];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("valid facet");

    let mut loads = [0.0_f64; 12];
    boundary
        .add_consistent_loads(&nodes, &mut loads)
        .expect("well-shaped");

    let total: f64 = (0..4).map(|node| loads[node * 3 + 1]).sum();
    let expected = traction[1] * 6.0;
    assert!((total - expected).abs() <= expected.abs() * f64::EPSILON * 8.0);
    // Node 3 is not on the facet and must carry nothing.
    for component in 0..3 {
        assert_eq!(loads[3 * 3 + component], 0.0);
    }
}

#[test]
fn a_traction_scales_with_the_facet_rather_than_the_node_count() {
    // The distinction between a traction and a force: doubling the facet
    // doubles the load. A load applied per node instead would not move.
    let short = [[0.0_f64, 0.0], [0.0, 1.0]];
    let long = [[0.0_f64, 0.0], [0.0, 2.0]];
    let traction = [3.0e5_f64, 0.0];
    let facets = [TractionFacet::new([0, 1], traction)];

    let mut short_loads = [0.0_f64; 4];
    TractionBoundary::try_new(&facets, &short)
        .expect("valid")
        .add_consistent_loads(&short, &mut short_loads)
        .expect("well-shaped");
    let mut long_loads = [0.0_f64; 4];
    TractionBoundary::try_new(&facets, &long)
        .expect("valid")
        .add_consistent_loads(&long, &mut long_loads)
        .expect("well-shaped");

    let short_total: f64 = (0..2).map(|node| short_loads[node * 2]).sum();
    let long_total: f64 = (0..2).map(|node| long_loads[node * 2]).sum();
    assert!((long_total - 2.0 * short_total).abs() <= short_total.abs() * f64::EPSILON * 8.0);
}

#[test]
fn tractions_accumulate_rather_than_overwrite() {
    // Two facets sharing a node, plus a pre-existing body force. Assignment
    // instead of accumulation would drop every load but the last.
    let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [2.0, 0.0]];
    let traction = [0.0_f64, 1.0e5];
    let facets = [
        TractionFacet::new([0, 1], traction),
        TractionFacet::new([1, 2], traction),
    ];
    let boundary = TractionBoundary::try_new(&facets, &nodes).expect("valid");

    let mut loads = [0.0_f64, 7.0, 0.0, 0.0, 0.0, 0.0];
    boundary
        .add_consistent_loads(&nodes, &mut loads)
        .expect("well-shaped");

    // Node 1 is shared, so it takes half of each facet: a full facet's worth.
    let expected_shared = traction[1] * 1.0;
    assert!((loads[3] - expected_shared).abs() <= expected_shared * f64::EPSILON * 8.0);
    // The pre-existing entry survived.
    assert!((loads[1] - (7.0 + traction[1] / 2.0)).abs() <= traction[1] * f64::EPSILON * 8.0);
}

#[test]
fn an_invalid_traction_boundary_is_rejected() {
    let nodes = [[0.0_f64, 0.0], [1.0, 0.0]];
    let facets = [TractionFacet::new([0, 5], [1.0_f64, 0.0])];
    assert_eq!(
        TractionBoundary::try_new(&facets, &nodes).expect_err("node 5 is beyond the mesh"),
        InvalidBoundary::NodeOutOfRange {
            facet: 0,
            position: 1,
            node: 5,
            nodes: 2,
        }
    );

    let facets = [TractionFacet::new([0, 0], [1.0_f64, 0.0])];
    assert_eq!(
        TractionBoundary::try_new(&facets, &nodes).expect_err("the facet has no length"),
        InvalidBoundary::DegenerateFacet { facet: 0 }
    );
}

// ---------------------------------------------------------------------------
// Dirichlet: validation
// ---------------------------------------------------------------------------

#[test]
fn dirichlet_conditions_are_validated_against_the_mesh() {
    let out_of_range = [PrescribedDisplacement::new(9, 0, 0.0_f64)];
    assert_eq!(
        DirichletConditions::<f64, 2>::try_new(&out_of_range, 5)
            .expect_err("node 9 is beyond the mesh"),
        InvalidConditions::NodeOutOfRange {
            position: 0,
            node: 9,
            nodes: 5,
        }
    );

    let bad_component = [PrescribedDisplacement::new(0, 2, 0.0_f64)];
    assert_eq!(
        DirichletConditions::<f64, 2>::try_new(&bad_component, 5)
            .expect_err("component 2 does not exist in 2-D"),
        InvalidConditions::ComponentOutOfRange {
            position: 0,
            component: 2,
            dimensions: 2,
        }
    );
}

#[test]
fn duplicate_and_unordered_conditions_are_rejected() {
    // A duplicated degree of freedom would resolve by whichever condition the
    // caller listed last — a silent dependence on input order.
    let duplicated = [
        PrescribedDisplacement::new(1, 0, 0.0_f64),
        PrescribedDisplacement::new(1, 0, 5.0),
    ];
    assert_eq!(
        DirichletConditions::<f64, 2>::try_new(&duplicated, 5)
            .expect_err("degree of freedom 2 appears twice"),
        InvalidConditions::NotStrictlyIncreasing {
            position: 1,
            degree_of_freedom: 2,
            previous: 2,
        }
    );

    let unordered = [
        PrescribedDisplacement::new(2, 0, 0.0_f64),
        PrescribedDisplacement::new(1, 0, 0.0),
    ];
    assert_eq!(
        DirichletConditions::<f64, 2>::try_new(&unordered, 5).expect_err("the conditions descend"),
        InvalidConditions::NotStrictlyIncreasing {
            position: 1,
            degree_of_freedom: 2,
            previous: 4,
        }
    );
}

// ---------------------------------------------------------------------------
// Dirichlet: the constrained operator
// ---------------------------------------------------------------------------

/// Fix both components of nodes 0 and 3, with node 3 prescribed non-zero.
fn conditions() -> [PrescribedDisplacement<f64>; 4] {
    [
        PrescribedDisplacement::new(0, 0, 0.0),
        PrescribedDisplacement::new(0, 1, 0.0),
        PrescribedDisplacement::new(3, 0, 2.5e-4),
        PrescribedDisplacement::new(3, 1, -1.1e-4),
    ]
}

#[test]
fn the_constrained_operator_is_symmetric() {
    // The defect this guards: projecting the output but not the input leaves
    // the coupling column intact, and the resulting operator still returns
    // plausible numbers while invalidating conjugate gradients.
    let (nodes, cells) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let prescribed = conditions();
    let bc = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");
    let m = moduli::<f64>(200e9, 0.3);

    let u: [f64; 10] = core::array::from_fn(|i| 1e-4 * ((i * 7 % 5) as f64 - 2.0));
    let v: [f64; 10] = core::array::from_fn(|i| 1e-4 * ((i * 3 % 7) as f64 - 3.0));
    let mut au = [0.0_f64; 10];
    let mut av = [0.0_f64; 10];
    let mut scratch = [0.0_f64; 10];
    mesh.constrained_action(&m, &bc, &u, &mut au, &mut scratch)
        .expect("well-shaped");
    mesh.constrained_action(&m, &bc, &v, &mut av, &mut scratch)
        .expect("well-shaped");

    let left: f64 = v.iter().zip(au.iter()).map(|(a, b)| a * b).sum();
    let right: f64 = u.iter().zip(av.iter()).map(|(a, b)| a * b).sum();
    let scale = left.abs().max(right.abs());
    assert!(
        (left - right).abs() <= scale * f64::EPSILON * 64.0,
        "the constrained operator is not symmetric: {left} vs {right}"
    );
}

#[test]
fn the_constrained_operator_is_positive_definite() {
    // Positive on *every* non-zero field, not merely non-rigid ones: the
    // conditions remove the rigid-body modes, which is what makes the
    // constrained system solvable at all. A rigid translation is now a
    // positive-energy field rather than a null one.
    let (nodes, cells) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let prescribed = conditions();
    let bc = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");
    let m = moduli::<f64>(200e9, 0.3);
    let mut scratch = [0.0_f64; 10];

    let translation = [1e-4_f64, 2e-4].repeat(5);
    let mut image = [0.0_f64; 10];
    mesh.constrained_action(&m, &bc, &translation, &mut image, &mut scratch)
        .expect("well-shaped");
    let energy: f64 = translation
        .iter()
        .zip(image.iter())
        .map(|(a, b)| a * b)
        .sum();
    assert!(
        energy > 0.0,
        "the constrained operator has a rigid translation in its null space: {energy}"
    );
}

#[test]
fn the_constrained_operator_leaves_identity_rows() {
    // The constrained part of the output is the constrained part of the input,
    // bit for bit. That is what makes the constrained rows read `u_c = g` once
    // the load carries `g`.
    let (nodes, cells) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let prescribed = conditions();
    let bc = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");
    let m = moduli::<f64>(200e9, 0.3);

    let u: [f64; 10] = core::array::from_fn(|i| 1e-4 * (i as f64 + 1.0));
    let mut image = [0.0_f64; 10];
    let mut scratch = [0.0_f64; 10];
    mesh.constrained_action(&m, &bc, &u, &mut image, &mut scratch)
        .expect("well-shaped");

    for condition in &prescribed {
        let dof = condition.degree_of_freedom::<2>();
        assert_eq!(image[dof], u[dof], "row {dof} is not an identity row");
    }
}

#[test]
fn the_constrained_operator_does_not_modify_its_input() {
    let (nodes, cells) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let prescribed = conditions();
    let bc = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");

    let u: [f64; 10] = core::array::from_fn(|i| 1e-4 * (i as f64 + 1.0));
    let before = u;
    let mut image = [0.0_f64; 10];
    let mut scratch = [0.0_f64; 10];
    mesh.constrained_action(
        &moduli::<f64>(200e9, 0.3),
        &bc,
        &u,
        &mut image,
        &mut scratch,
    )
    .expect("well-shaped");
    assert_eq!(u, before);
}

#[test]
fn the_constrained_load_carries_the_prescribed_displacement() {
    // The load's constrained rows must hold `g` itself, so that the identity
    // rows of the operator solve to `u_c = g`.
    let (nodes, cells) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let prescribed = conditions();
    let bc = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");

    let external = [0.0_f64; 10];
    let mut load = [0.0_f64; 10];
    let mut scratch = [0.0_f64; 10];
    mesh.constrained_load(
        &moduli::<f64>(200e9, 0.3),
        &bc,
        &external,
        &mut load,
        &mut scratch,
    )
    .expect("well-shaped");

    for condition in &prescribed {
        assert_eq!(
            load[condition.degree_of_freedom::<2>()],
            *condition.value(),
            "the load does not carry the prescribed value"
        );
    }
}

#[test]
fn a_non_zero_prescribed_displacement_reaches_the_free_rows() {
    // The defect this guards is the one that survives most test suites:
    // dropping `K u_g` from the load silently solves the problem in which
    // every prescribed displacement is zero. With the conditions all zero the
    // two loads agree, so only a non-zero prescription separates them.
    let (nodes, cells) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let m = moduli::<f64>(200e9, 0.3);
    let external = [0.0_f64; 10];
    let mut scratch = [0.0_f64; 10];

    let moved = conditions();
    let bc = DirichletConditions::try_new(&moved, mesh.node_count()).expect("valid");
    let mut load = [0.0_f64; 10];
    mesh.constrained_load(&m, &bc, &external, &mut load, &mut scratch)
        .expect("well-shaped");

    let held = [
        PrescribedDisplacement::new(0, 0, 0.0_f64),
        PrescribedDisplacement::new(0, 1, 0.0),
        PrescribedDisplacement::new(3, 0, 0.0),
        PrescribedDisplacement::new(3, 1, 0.0),
    ];
    let zeroed = DirichletConditions::try_new(&held, mesh.node_count()).expect("valid");
    let mut zero_load = [0.0_f64; 10];
    mesh.constrained_load(&m, &zeroed, &external, &mut zero_load, &mut scratch)
        .expect("well-shaped");

    // The free rows must differ: the prescribed motion strains the cells that
    // touch node 3, and that strain is a load on their other nodes.
    let free: Vec<usize> = (0..10)
        .filter(|dof| !moved.iter().any(|c| c.degree_of_freedom::<2>() == *dof))
        .collect();
    let difference: f64 = free
        .iter()
        .map(|dof| (load[*dof] - zero_load[*dof]).abs())
        .fold(0.0, f64::max);
    assert!(
        difference > 0.0,
        "the prescribed displacement did not reach any free row, so the load is solving the \
         all-zero problem"
    );

    // And the zero case must itself be zero, confirming the difference above
    // comes from the prescription rather than from noise in both.
    for dof in &free {
        assert_eq!(zero_load[*dof], 0.0);
    }
}

#[test]
fn misshaped_buffers_are_rejected_by_the_constrained_operator() {
    let (nodes, cells) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let prescribed = conditions();
    let bc = DirichletConditions::try_new(&prescribed, mesh.node_count()).expect("valid");
    let m = moduli::<f64>(200e9, 0.3);

    let mut image = [0.0_f64; 10];
    let mut short = [0.0_f64; 4];
    assert!(
        mesh.constrained_action(&m, &bc, &[0.0; 10], &mut image, &mut short)
            .is_err()
    );
    let mut scratch = [0.0_f64; 10];
    assert!(
        mesh.constrained_action(&m, &bc, &[0.0; 6], &mut image, &mut scratch)
            .is_err()
    );
}
