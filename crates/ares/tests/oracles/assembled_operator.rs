//! Properties of the assembled operator that the Krylov solve depends on.

use super::support::{moduli, square_patch};
use ares::SimplexMesh;

// ---------------------------------------------------------------------------
// Operator properties the solver depends on
// ---------------------------------------------------------------------------

#[test]
fn the_assembled_operator_is_symmetric() {
    // v . K u == u . K v. Conjugate gradients is only valid on a symmetric
    // operator, so an asymmetric assembly would not fail loudly — it would
    // converge to the wrong displacement or stall.
    let (nodes, cells, _) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let m = moduli::<f64>(200e9, 0.3);

    let u: [f64; 10] = core::array::from_fn(|i| 1e-4 * ((i * 7 % 5) as f64 - 2.0));
    let v: [f64; 10] = core::array::from_fn(|i| 1e-4 * ((i * 3 % 7) as f64 - 3.0));
    let mut ku = [0.0_f64; 10];
    let mut kv = [0.0_f64; 10];
    mesh.internal_forces(&m, &u, &mut ku).expect("well-shaped");
    mesh.internal_forces(&m, &v, &mut kv).expect("well-shaped");

    let left: f64 = v.iter().zip(ku.iter()).map(|(a, b)| a * b).sum();
    let right: f64 = u.iter().zip(kv.iter()).map(|(a, b)| a * b).sum();
    let scale = left.abs().max(right.abs());
    assert!(
        (left - right).abs() <= scale * f64::EPSILON * 64.0,
        "the operator is not symmetric: {left} vs {right}"
    );
}

#[test]
fn the_assembled_operator_is_positive_on_non_rigid_modes() {
    // u . K u > 0 away from the rigid-body null space. Conjugate gradients
    // requires positive definiteness on the constrained subspace, and a single
    // inverted cell would break it.
    let (nodes, cells, _) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let m = moduli::<f64>(200e9, 0.3);

    let u: [f64; 10] = core::array::from_fn(|i| 1e-4 * ((i % 3) as f64 - 1.0) * (i as f64 + 1.0));
    let mut ku = [0.0_f64; 10];
    mesh.internal_forces(&m, &u, &mut ku).expect("well-shaped");
    let energy: f64 = u.iter().zip(ku.iter()).map(|(a, b)| a * b).sum();
    assert!(
        energy > 0.0,
        "a strained mesh stored {energy}, expected > 0"
    );
}

#[test]
fn the_assembled_operator_is_linear() {
    let (nodes, cells, _) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let m = moduli::<f64>(200e9, 0.3);
    let base: [f64; 10] = core::array::from_fn(|i| 1e-4 * ((i % 4) as f64 - 1.5));
    let factor = -2.75_f64;
    let scaled: [f64; 10] = core::array::from_fn(|i| base[i] * factor);

    let mut single = [0.0_f64; 10];
    let mut multiple = [0.0_f64; 10];
    mesh.internal_forces(&m, &base, &mut single)
        .expect("well-shaped");
    mesh.internal_forces(&m, &scaled, &mut multiple)
        .expect("well-shaped");
    for (dof, (one, many)) in single.iter().zip(multiple.iter()).enumerate() {
        let expected = one * factor;
        let tolerance = expected.abs().max(1.0) * f64::EPSILON * 32.0;
        assert!(
            (many - expected).abs() <= tolerance,
            "dof {dof}: {many} != {expected}"
        );
    }
}
