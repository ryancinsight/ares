//! Element geometry: measures, shape gradients, and their identities.
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

use super::support::{moduli, unit_tetrahedron, unit_triangle};
use ares::{DegenerateElement, Simplex, stiffness_action};

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
