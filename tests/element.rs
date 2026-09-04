//! Executable evidence for the Phase 0 element kernel (atlas ADR 0057).
//!
//! The load-bearing oracle here is the rigid-body null space: an element must
//! exert no force under a rigid motion. It is the element-level precursor to
//! the patch test, and the failure it guards is silent — an element with a
//! leaky null space stiffens rigid motion, producing a spurious stress that
//! grows with how far the body has moved rather than with the mesh, so
//! refinement never reveals it.
//!
//! The two rigid motions are asserted differently, and the difference is real
//! rather than a matter of care. **Translation is exact**, because
//! `stiffness_action` differences displacements against a reference node, so
//! every relative displacement is identically zero. **Rotation is bounded**,
//! because its relative displacements are not zero and the reconstructed
//! gradient is antisymmetric only to rounding. Asserting exactness for
//! rotation would assert a coincidence of one geometry: it does hold for the
//! unit triangle, and does not in general.

#![expect(
    clippy::float_cmp,
    reason = "the exact comparisons are exact by construction: unit-simplex measures are ratios of small integers, and translation-invariant forces are identically zero because the relative displacements feeding them are. Everything with genuine rounding - the partition of unity, the linear-field reproduction, rotation forces, and equilibrium - carries a derived bound instead."
)]

use aequitas::systems::si::quantities::{Dimensionless, Pressure};
use ares::{DegenerateElement, Simplex, stiffness_action};
use eunomia::RealField;
use proteus::IsotropicModuli;

fn moduli<T: RealField>(young: f64, poisson: f64) -> IsotropicModuli<T> {
    IsotropicModuli::from_young_poisson(
        Pressure::from_base(T::from_f64(young)),
        Dimensionless::from_base(T::from_f64(poisson)),
    )
    .expect("inside the positive-definite domain")
}

/// The unit triangle: (0,0), (1,0), (0,1).
fn unit_triangle() -> [[f64; 2]; 3] {
    [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
}

/// The unit tetrahedron: origin plus the three axis points.
fn unit_tetrahedron() -> [[f64; 3]; 4] {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

#[test]
fn unit_simplex_measures_match_their_closed_forms() {
    let nodes = unit_triangle();
    let triangle = Simplex::new(&nodes);
    assert_eq!(triangle.signed_measure(), 0.5);

    let nodes = unit_tetrahedron();
    let tetrahedron = Simplex::new(&nodes);
    assert_eq!(tetrahedron.signed_measure(), 1.0 / 6.0);
}

#[test]
fn reversing_node_order_flips_the_sign_of_the_measure() {
    // A negative measure means an inverted element, not a small one. Silently
    // taking the absolute value would hide an inverted mesh.
    let nodes = [[0.0_f64, 0.0], [0.0, 1.0], [1.0, 0.0]];
    let flipped = Simplex::new(&nodes);
    assert_eq!(flipped.signed_measure(), -0.5);
}

#[test]
fn shape_gradients_cancel_to_rounding() {
    // Partition of unity. `grad N_0` is the negated sum of the others, so the
    // cancellation is exact only when re-summed in that same order; summed in
    // another it leaves a residual, measured at 1.4e-17 for an ordinary
    // triangle. An earlier version of this test asserted exact zero and passed
    // — on a geometry that happens to cancel. The bound is what actually
    // holds, and `stiffness_action` does not depend on the stronger claim.
    let nodes = [[0.3_f64, -1.2], [2.7, 0.4], [-0.9, 3.1]];
    let element = Simplex::new(&nodes);
    let gradients = element.shape_gradients().expect("non-degenerate");

    for component in 0..2 {
        let total: f64 = gradients.iter().map(|g| g[component]).sum();
        let scale = gradients
            .iter()
            .map(|g| g[component].abs())
            .fold(1.0_f64, f64::max);
        assert!(
            total.abs() <= scale * f64::EPSILON * 8.0,
            "component {component} residual {total} exceeds the rounding bound"
        );
    }
}

#[test]
fn shape_gradients_reproduce_a_linear_field() {
    // Consistency: for a linear field f(x) = a.x, sum_i f(x_i) grad N_i = a.
    // This is the property that makes constant strain representable, which the
    // patch test then exercises globally.
    let nodes = [[0.2_f64, 0.1], [1.7, 0.3], [0.4, 2.2]];
    let element = Simplex::new(&nodes);
    let gradients = element.shape_gradients().expect("non-degenerate");

    let a = [1.7_f64, -0.6];
    let mut recovered = [0.0_f64; 2];
    for (node, gradient) in nodes.iter().zip(gradients.iter()) {
        let value = a[0] * node[0] + a[1] * node[1];
        for (component, slot) in recovered.iter_mut().enumerate() {
            *slot += value * gradient[component];
        }
    }
    for (component, slot) in recovered.iter().enumerate() {
        let tolerance = a[component].abs() * f64::EPSILON * 32.0;
        assert!(
            (*slot - a[component]).abs() <= tolerance,
            "component {component}: {slot} != {}",
            a[component]
        );
    }
}

#[test]
fn a_degenerate_element_is_rejected_rather_than_returning_infinities() {
    // Three collinear nodes have no area, so no shape gradient exists.
    // Returning infinities would push a plausible NaN into the assembled
    // system, where it is far harder to attribute.
    let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [2.0, 0.0]];
    let element = Simplex::new(&nodes);
    assert_eq!(element.shape_gradients(), Err(DegenerateElement));
}

// ---------------------------------------------------------------------------
// Rigid-body null space
// ---------------------------------------------------------------------------

fn assert_translation_is_force_free<T: RealField>() {
    let nodes = [
        [T::from_f64(0.1), T::from_f64(0.2)],
        [T::from_f64(1.3), T::from_f64(0.1)],
        [T::from_f64(0.4), T::from_f64(1.9)],
    ];
    let element = Simplex::new(&nodes);
    let shift = [T::from_f64(3.7), T::from_f64(-2.1)];
    let displacements = [shift, shift, shift];
    let mut forces = [[T::from_f64(0.0); 2]; 3];

    stiffness_action(
        &element,
        &moduli::<T>(200e9, 0.3),
        &displacements,
        &mut forces,
    )
    .expect("valid element");

    for force in &forces {
        for component in force {
            assert!(
                *component == T::from_f64(0.0),
                "translation produced a nodal force"
            );
        }
    }
}

#[test]
fn rigid_translation_produces_exactly_zero_nodal_forces() {
    assert_translation_is_force_free::<f32>();
    assert_translation_is_force_free::<f64>();
}

#[test]
fn infinitesimal_rotation_of_the_unit_triangle_is_force_free() {
    // Exact for *this* geometry, whose edge matrix is the identity so the
    // shape gradients carry no rounding. The general case is bounded rather
    // than exact — see the rotation property test.
    let nodes = unit_triangle();
    let element = Simplex::new(&nodes);
    let w = 0.42_f64;
    let displacements = [
        [-w * nodes[0][1], w * nodes[0][0]],
        [-w * nodes[1][1], w * nodes[1][0]],
        [-w * nodes[2][1], w * nodes[2][0]],
    ];
    let mut forces = [[0.0_f64; 2]; 3];
    stiffness_action(
        &element,
        &moduli::<f64>(70e9, 0.33),
        &displacements,
        &mut forces,
    )
    .expect("valid element");

    for force in &forces {
        for component in force {
            assert_eq!(*component, 0.0, "rotation produced a nodal force");
        }
    }
}

#[test]
fn a_rigid_motion_of_a_tetrahedron_is_also_force_free() {
    let nodes = unit_tetrahedron();
    let element = Simplex::new(&nodes);
    let shift = [1.5_f64, -0.25, 3.0];
    let displacements = [shift, shift, shift, shift];
    let mut forces = [[0.0_f64; 3]; 4];
    stiffness_action(
        &element,
        &moduli::<f64>(200e9, 0.3),
        &displacements,
        &mut forces,
    )
    .expect("valid element");

    for force in &forces {
        for component in force {
            assert_eq!(*component, 0.0);
        }
    }
}

// ---------------------------------------------------------------------------
// Equilibrium and energy
// ---------------------------------------------------------------------------

#[test]
fn nodal_forces_are_self_equilibrated() {
    // Newton's third law at element level: with no body force, the nodal
    // forces an element exerts must sum to zero for any displacement, rigid or
    // not. This follows from the shape gradients summing to zero, so it is a
    // second consumer of that identity.
    let nodes = [[0.0_f64, 0.0], [2.0, 0.3], [0.5, 1.7]];
    let element = Simplex::new(&nodes);
    let displacements = [[0.01_f64, -0.02], [0.03, 0.005], [-0.015, 0.02]];
    let mut forces = [[0.0_f64; 2]; 3];
    stiffness_action(
        &element,
        &moduli::<f64>(200e9, 0.3),
        &displacements,
        &mut forces,
    )
    .expect("valid element");

    for component in 0..2 {
        let total: f64 = forces.iter().map(|f| f[component]).sum();
        let scale = forces
            .iter()
            .map(|f| f[component].abs())
            .fold(1.0_f64, f64::max);
        assert!(
            total.abs() <= scale * f64::EPSILON * 16.0,
            "component {component} is not equilibrated: {total}"
        );
    }
}

#[test]
fn the_stiffness_action_is_energetically_positive() {
    // u . K u > 0 for a non-rigid displacement: the element stores energy
    // under deformation. A negative value would mean an element that releases
    // energy when strained, which makes the global solve indefinite.
    let nodes = unit_triangle();
    let element = Simplex::new(&nodes);
    let displacements = [[0.0_f64, 0.0], [0.01, 0.0], [0.0, 0.0]];
    let mut forces = [[0.0_f64; 2]; 3];
    stiffness_action(
        &element,
        &moduli::<f64>(200e9, 0.3),
        &displacements,
        &mut forces,
    )
    .expect("valid element");

    let energy: f64 = displacements
        .iter()
        .zip(forces.iter())
        .map(|(u, f)| u[0] * f[0] + u[1] * f[1])
        .sum();
    assert!(
        energy > 0.0,
        "strained element stored {energy}, expected > 0"
    );
}

#[test]
fn the_action_is_linear_in_displacement() {
    // K(k u) = k K(u): linearity is what makes the assembled operator
    // constant, so a nonlinearity here would silently invalidate every solve.
    let nodes = unit_triangle();
    let element = Simplex::new(&nodes);
    let m = moduli::<f64>(200e9, 0.3);
    let base = [[0.0_f64, 0.0], [0.002, 0.001], [-0.001, 0.003]];
    let k = 3.5_f64;
    let scaled = base.map(|u| [u[0] * k, u[1] * k]);

    let mut single = [[0.0_f64; 2]; 3];
    let mut multiple = [[0.0_f64; 2]; 3];
    stiffness_action(&element, &m, &base, &mut single).expect("valid");
    stiffness_action(&element, &m, &scaled, &mut multiple).expect("valid");

    for (one, many) in single.iter().zip(multiple.iter()) {
        for component in 0..2 {
            let expected = one[component] * k;
            let tolerance = expected.abs().max(1.0) * f64::EPSILON * 16.0;
            assert!((many[component] - expected).abs() <= tolerance);
        }
    }
}

// ---------------------------------------------------------------------------
// Properties
// ---------------------------------------------------------------------------

proptest::proptest! {
    #[test]
    fn any_rigid_translation_is_force_free(
        ux in -1e3_f64..1e3,
        uy in -1e3_f64..1e3,
        young in 1e6_f64..5e11,
    ) {
        let nodes = [[0.0_f64, 0.0], [1.0, 0.0], [0.0, 1.0]];
        let element = Simplex::new(&nodes);
        let shift = [ux, uy];
        let displacements = [shift, shift, shift];
        let mut forces = [[0.0_f64; 2]; 3];
        stiffness_action(&element, &moduli::<f64>(young, 0.3), &displacements, &mut forces)
            .expect("valid element");

        for force in &forces {
            for component in force {
                proptest::prop_assert_eq!(*component, 0.0);
            }
        }
    }

    #[test]
    fn any_rigid_rotation_leaves_forces_at_rounding(
        w in -1.0_f64..1.0,
        x in 0.5_f64..3.0,
        y in 0.5_f64..3.0,
        young in 1e6_f64..5e11,
    ) {
        // Rotation is *not* exact the way translation is. Translation is exact
        // by construction because every (u_a - u_0) vanishes; a rotation's
        // differences do not, so the reconstructed gradient is antisymmetric
        // only to rounding and the strain it feeds is small rather than zero.
        // Asserting exactness here would be asserting a coincidence of one
        // geometry — this bounds it instead.
        let nodes = [[0.0_f64, 0.0], [x, 0.0], [0.0, y]];
        let element = Simplex::new(&nodes);
        let displacements = nodes.map(|p| [-w * p[1], w * p[0]]);
        let mut forces = [[0.0_f64; 2]; 3];
        stiffness_action(&element, &moduli::<f64>(young, 0.3), &displacements, &mut forces)
            .expect("valid element");

        // Scale: stiffness times the largest displacement in the element.
        let span = x.max(y) * w.abs();
        let bound = young * span * f64::EPSILON * 64.0;
        for force in &forces {
            for component in force {
                proptest::prop_assert!(
                    component.abs() <= bound,
                    "rotation force {} exceeds the rounding bound {bound}", component
                );
            }
        }
    }

    #[test]
    fn shape_gradients_always_cancel_to_rounding(
        x in 0.5_f64..3.0, y in 0.5_f64..3.0, skew in -2.0_f64..2.0,
    ) {
        // Skewed rather than axis-aligned: an axis-aligned triangle has a
        // diagonal edge matrix and cancels exactly, which would make this
        // property test agree with the stronger claim for the wrong reason.
        let nodes = [[0.0_f64, 0.0], [x, skew], [skew, y]];
        let element = Simplex::new(&nodes);
        let Ok(gradients) = element.shape_gradients() else {
            return Ok(());
        };
        for component in 0..2 {
            let total: f64 = gradients.iter().map(|g| g[component]).sum();
            let scale = gradients
                .iter()
                .map(|g| g[component].abs())
                .fold(1.0_f64, f64::max);
            proptest::prop_assert!(total.abs() <= scale * f64::EPSILON * 8.0);
        }
    }
}

// ---------------------------------------------------------------------------
// Hand-computed element stiffness
// ---------------------------------------------------------------------------

/// Lame parameters from the engineering pair.
fn lame(young: f64, poisson: f64) -> (f64, f64) {
    (
        young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson)),
        young / (2.0 * (1.0 + poisson)),
    )
}

/// The stiffness column the action produces for a unit displacement at
/// `(node, component)`.
fn stiffness_column<const D: usize, const N: usize>(
    nodes: &[[f64; D]; N],
    young: f64,
    poisson: f64,
    node: usize,
    component: usize,
) -> [[f64; D]; N] {
    let mut displacements = [[0.0_f64; D]; N];
    displacements[node][component] = 1.0;
    let mut forces = [[0.0_f64; D]; N];
    stiffness_action(
        &Simplex::new(nodes),
        &moduli::<f64>(young, poisson),
        &displacements,
        &mut forces,
    )
    .expect("valid element");
    forces
}

fn assert_column(actual: &[[f64; 2]; 3], expected: &[[f64; 2]; 3], label: &str) {
    for (node, (got, want)) in actual.iter().zip(expected.iter()).enumerate() {
        for component in 0..2 {
            let tolerance = want[component].abs().max(1.0) * f64::EPSILON * 16.0;
            assert!(
                (got[component] - want[component]).abs() <= tolerance,
                "{label}: node {node} component {component} is {} but hand computation gives {}",
                got[component],
                want[component]
            );
        }
    }
}

#[test]
fn the_unit_triangle_matches_its_hand_computed_stiffness_columns() {
    // Derived independently of the implementation. The reference route is the
    // textbook Voigt one — `K = A B^T D B` with
    //
    //   B rows  [-1, 0, 1, 0, 0, 0]        (eps_xx)
    //           [ 0,-1, 0, 0, 0, 1]        (eps_yy)
    //           [-1,-1, 0, 1, 1, 0]        (gamma_xy, the engineering shear)
    //   D       [[l+2m, l, 0], [l, l+2m, 0], [0, 0, m]]
    //   A       1/2
    //
    // while the implementation goes through the full tensor and never forms a
    // `B`. Agreement across those two routes is what makes this an oracle
    // rather than a restatement: the engineering-shear factor of two that
    // separates `gamma_xy` from `eps_xy` is present in the reference and
    // absent from the implementation, so a version that confused them would
    // disagree here.
    let (young, poisson) = (200e9, 0.3);
    let (l, m) = lame(young, poisson);
    let nodes = unit_triangle();
    let half = 0.5;

    // Column for u_0x: A * [l+3m, l+m, -(l+2m), -m, -m, -l].
    assert_column(
        &stiffness_column(&nodes, young, poisson, 0, 0),
        &[
            [half * (l + 3.0 * m), half * (l + m)],
            [half * -(l + 2.0 * m), half * -m],
            [half * -m, half * -l],
        ],
        "column u_0x",
    );

    // Column for u_0y: A * [l+m, l+3m, -l, -m, -m, -(l+2m)].
    assert_column(
        &stiffness_column(&nodes, young, poisson, 0, 1),
        &[
            [half * (l + m), half * (l + 3.0 * m)],
            [half * -l, half * -m],
            [half * -m, half * -(l + 2.0 * m)],
        ],
        "column u_0y",
    );

    // Column for u_1x: A * [-(l+2m), -l, l+2m, 0, 0, l].
    assert_column(
        &stiffness_column(&nodes, young, poisson, 1, 0),
        &[
            [half * -(l + 2.0 * m), half * -l],
            [half * (l + 2.0 * m), 0.0],
            [0.0, half * l],
        ],
        "column u_1x",
    );
}

#[test]
fn the_unit_tetrahedron_matches_its_hand_computed_stiffness_column() {
    // The same derivation in 3-D for u_1x. With grad N_1 = (1,0,0) the
    // displacement gradient is a single unit entry, so
    //   eps = diag(1,0,0), tr = 1, sigma = diag(l+2m, l, l)
    // and f_a = V sigma . grad N_a with V = 1/6.
    let (young, poisson) = (200e9, 0.3);
    let (l, m) = lame(young, poisson);
    let nodes = unit_tetrahedron();
    let volume = 1.0 / 6.0;

    let actual = stiffness_column(&nodes, young, poisson, 1, 0);
    let expected = [
        [volume * -(l + 2.0 * m), volume * -l, volume * -l],
        [volume * (l + 2.0 * m), 0.0, 0.0],
        [0.0, volume * l, 0.0],
        [0.0, 0.0, volume * l],
    ];
    for (node, (got, want)) in actual.iter().zip(expected.iter()).enumerate() {
        for component in 0..3 {
            let tolerance = want[component].abs().max(1.0) * f64::EPSILON * 16.0;
            assert!(
                (got[component] - want[component]).abs() <= tolerance,
                "node {node} component {component} is {} but hand computation gives {}",
                got[component],
                want[component]
            );
        }
    }
}

#[test]
fn the_element_stiffness_is_symmetric() {
    // K = K^T across every column pair, which the three hand-computed columns
    // above check only where they overlap.
    let (young, poisson) = (200e9, 0.3);
    let nodes = [[0.2_f64, 0.1], [1.7, 0.3], [0.4, 2.2]];
    let mut matrix = [[0.0_f64; 6]; 6];
    for node in 0..3 {
        for component in 0..2 {
            let column = stiffness_column(&nodes, young, poisson, node, component);
            for (row_node, force) in column.iter().enumerate() {
                for (row_component, value) in force.iter().enumerate() {
                    matrix[row_node * 2 + row_component][node * 2 + component] = *value;
                }
            }
        }
    }
    let scale = matrix
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    for (i, row) in matrix.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            let transposed = matrix[j][i];
            assert!(
                (value - transposed).abs() <= scale * f64::EPSILON * 16.0,
                "K[{i}][{j}] = {value} but K[{j}][{i}] = {transposed}"
            );
        }
    }
}
