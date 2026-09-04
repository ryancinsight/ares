//! Dirichlet conditions and the constrained operator.
//!
//! The constrained operator must stay symmetric - conjugate gradients has no
//! convergence guarantee otherwise, and the failure is a wrong answer rather
//! than an error. The two defects guarded here are both silent: projecting on
//! one side only leaves `K_cf` intact and breaks symmetry, and dropping the
//! `K u_g` term from the load silently solves the problem in which every
//! prescribed displacement is zero. That second one is why a non-zero
//! prescribed value appears in these tests at all - fixed-at-zero is both the
//! common case and the one where the omission is undetectable.

use super::support::{moduli, square_patch};
use ares::{DirichletConditions, InvalidConditions, PrescribedDisplacement, SimplexMesh};

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
    let (nodes, cells, _) = square_patch();
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
    let (nodes, cells, _) = square_patch();
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
    let (nodes, cells, _) = square_patch();
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
    let (nodes, cells, _) = square_patch();
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
    let (nodes, cells, _) = square_patch();
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
    let (nodes, cells, _) = square_patch();
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
    let (nodes, cells, _) = square_patch();
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
