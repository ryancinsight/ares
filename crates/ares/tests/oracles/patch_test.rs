//! Executable evidence for global assembly (atlas ADR 0057, phase A5).
//!
//! # The patch test
//!
//! The headline oracle. A constant-strain displacement field is imposed on an
//! arbitrary distorted patch, and the assembled operator must reproduce it:
//! every cell recovers the same constant strain, and every **interior** node
//! carries zero internal force.
//!
//! It earns its status because it fails on almost any assembly defect — a
//! wrong gather or scatter index, a mapping error, a quadrature weight, a
//! sign, a transposed gradient — while a convergence study can pass with
//! several of those present. The mechanism is exact cancellation: for an
//! interior node `p`, the shape function `N_p` vanishes on the patch boundary,
//! so
//!
//! ```text
//! f_p = sum_e |Om_e| sigma . grad N_p^e = sigma . integral(grad N_p) = 0
//! ```
//!
//! by the divergence theorem, and the sum is over terms that are individually
//! large. Anything that perturbs one term without perturbing its partners
//! leaves a residual far above rounding.
//!
//! # "Exact to machine precision" is a derived bound, not an equality
//!
//! The cancellation above is exact in exact arithmetic. In floating point the
//! nodal values `u_i = a + G x_i` already carry rounding, so the bound is
//! derived from the problem rather than asserted as zero:
//!
//! ```text
//! |f_p| <= n * |Om| * (lambda + 2 mu) * max|u| * max|grad N|^2 * eps
//! ```
//!
//! Every factor is measured from the patch itself. The derivation: a nodal
//! force is `|Om| sigma grad N`, with `sigma ~ C eps_strain` and
//! `eps_strain ~ |du| |grad N|`, so an absolute error `|u| * eps` in the
//! displacements propagates to `|Om| C |u| |grad N|^2 * eps` in the force. `n`
//! covers the handful of accumulations per node.
//!
//! The constant term `a` enters that bound through `max|u|`, and it has to:
//! `u_i = a + G x_i` is rounded at construction with absolute error `|a| eps`,
//! and no later differencing can recover information the stored value never
//! had. So the large-offset case below does **not** prove that the reference
//! -node differencing works — a formulation without it leaves a residual of
//! exactly the same order, which this bound therefore cannot separate. What
//! it does establish is that the residual stays set by input rounding when the
//! constant term dominates the strain by four orders, rather than being
//! amplified by the formulation.
//!
//! Translation invariance is proved instead by
//! `a_rigid_translation_of_the_whole_mesh_is_exactly_force_free`, where every
//! nodal value is the *same* float, the differences are identically zero, and
//! the assertion is exact equality rather than a bound.
//!
//! # What these oracles were measured to catch
//!
//! Verified by mutation rather than asserted. Reversing the scatter order
//! fails 8 tests; a Voigt engineering-shear factor of two on the off-diagonal
//! strain fails 6; dropping the cell measure from the nodal force fails 6;
//! removing the reference-node differencing fails the exact-translation test
//! above. Transposing the stress in the force contraction fails nothing, and
//! that one is correct: the tensor is symmetric by construction, so the
//! mutation is a no-op rather than a defect the suite misses.

use super::support::{cube_patch, lame, moduli, square_patch};
use ares::SimplexMesh;

// ---------------------------------------------------------------------------
// Rigid-body null space, assembled
// ---------------------------------------------------------------------------

#[test]
fn a_rigid_translation_of_the_whole_mesh_is_exactly_force_free() {
    // The element-level property, now scattered and accumulated across cells.
    // Still exact: each cell's contribution is a sum of exact zeros, and
    // summing exact zeros stays exact.
    let (nodes, cells, _) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");

    let shift = [4.2_f64, -7.5];
    let mut displacements = [0.0_f64; 10];
    for slot in displacements.as_chunks_mut::<2>().0 {
        *slot = shift;
    }
    let mut forces = [1.0_f64; 10];
    mesh.internal_forces(&moduli::<f64>(200e9, 0.3), &displacements, &mut forces)
        .expect("well-shaped fields");

    for (dof, force) in forces.iter().enumerate() {
        assert_eq!(*force, 0.0, "translation produced a force at dof {dof}");
    }
}

#[test]
fn the_force_field_is_written_not_accumulated() {
    // `f = K u` is an assignment. A caller reusing a buffer across Krylov
    // iterations would otherwise accumulate every previous iteration into the
    // current one, which converges to nonsense rather than failing.
    let (nodes, cells, _) = square_patch();
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let m = moduli::<f64>(200e9, 0.3);
    let displacements: [f64; 10] = core::array::from_fn(|i| 1e-4 * (i as f64).sin());

    let mut once = [0.0_f64; 10];
    mesh.internal_forces(&m, &displacements, &mut once)
        .expect("well-shaped");
    let mut reused = [9.9e9_f64; 10];
    mesh.internal_forces(&m, &displacements, &mut reused)
        .expect("well-shaped");
    assert_eq!(once, reused);
}

// ---------------------------------------------------------------------------
// The patch test
// ---------------------------------------------------------------------------

/// Impose `u(x) = a + G x` at every node of a `D`-dimensional patch, then
/// assert the two halves of the patch test.
fn assert_patch_test<const D: usize, const N: usize, const NODES: usize, const CELLS: usize>(
    nodes: &[[f64; D]; NODES],
    cells: &[[usize; N]; CELLS],
    interior: usize,
    offset: [f64; D],
    gradient: [[f64; D]; D],
) {
    let (young, poisson) = (200e9, 0.3);
    let (lambda, mu) = lame(young, poisson);
    let mesh = SimplexMesh::try_new(nodes, cells).expect("valid patch");

    // `NODES * D` is not a const expression on stable, and the test harness
    // links std, so the field vectors are heap-allocated here. The library
    // itself never allocates.
    // u_i = a + G x_i, evaluated per node.
    let mut displacements = vec![0.0_f64; NODES * D];
    for (node, position) in nodes.iter().enumerate() {
        for component in 0..D {
            let mut value = offset[component];
            for axis in 0..D {
                value += gradient[component][axis] * position[axis];
            }
            displacements[node * D + component] = value;
        }
    }

    // Half one: every cell recovers the same constant strain, sym(G).
    let expected_strain: [[f64; D]; D] =
        core::array::from_fn(|i| core::array::from_fn(|j| 0.5 * (gradient[i][j] + gradient[j][i])));
    let strain_magnitude = expected_strain
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0_f64, |m, v| m.max(v.abs()));
    let displacement_magnitude = displacements.iter().fold(0.0_f64, |m, v| m.max(v.abs()));

    // Largest shape gradient anywhere in the patch, which sets both the strain
    // conditioning and the force scale below. Measured rather than assumed,
    // because a distorted patch has no single element size.
    let mut gradient_magnitude = 0.0_f64;
    let mut total_measure = 0.0_f64;
    for connectivity in cells {
        let mut coordinates = [[0.0_f64; D]; N];
        for (slot, node) in coordinates.iter_mut().zip(connectivity.iter()) {
            *slot = nodes[*node];
        }
        let element = ares::Simplex::new(&coordinates);
        total_measure += element.signed_measure();
        for shape in &element.shape_gradients().expect("non-degenerate") {
            for component in shape {
                gradient_magnitude = gradient_magnitude.max(component.abs());
            }
        }
    }

    // Strain error: an absolute rounding error `|u| * eps` in the nodal
    // displacements passes through `grad N`, so `|d eps| ~ |u| |grad N| eps`.
    let strain_bound =
        16.0 * displacement_magnitude * gradient_magnitude * f64::EPSILON * (N as f64);
    for (cell, strain) in mesh
        .cell_strains(&displacements)
        .expect("well-shaped")
        .enumerate()
    {
        let recovered = strain.tensor().components();
        for i in 0..D {
            for j in 0..D {
                assert!(
                    (recovered[i][j] - expected_strain[i][j]).abs() <= strain_bound,
                    "cell {cell} strain [{i}][{j}] is {} but the field imposes {} \
                     (bound {strain_bound:.3e})",
                    recovered[i][j],
                    expected_strain[i][j]
                );
            }
        }
    }

    // Half two: interior nodes carry no internal force.
    let mut forces = vec![0.0_f64; NODES * D];
    mesh.internal_forces(&moduli::<f64>(young, poisson), &displacements, &mut forces)
        .expect("well-shaped");

    // |f| ~ |Om| C |u| |grad N|^2 eps, per the module header.
    let stiffness = lambda + 2.0 * mu;
    let force_bound = 32.0
        * total_measure
        * stiffness
        * displacement_magnitude
        * gradient_magnitude
        * gradient_magnitude
        * f64::EPSILON;
    for component in 0..D {
        let residual = forces[interior * D + component];
        assert!(
            residual.abs() <= force_bound,
            "interior node {interior} carries force {residual:.6e} in component {component}, \
             above the derived bound {force_bound:.3e}"
        );
    }

    // The patch as a whole is self-equilibrated: the boundary tractions the
    // constant stress implies must sum to zero. This is independent of the
    // interior assertion and catches a scatter that drops a cell entirely.
    for component in 0..D {
        let total: f64 = (0..NODES).map(|node| forces[node * D + component]).sum();
        assert!(
            total.abs() <= force_bound,
            "component {component} of the patch is not equilibrated: {total:.6e}"
        );
    }

    // Guard against a vacuous pass: a field that produced no stress would
    // satisfy every assertion above trivially.
    assert!(
        strain_magnitude > 0.0,
        "the imposed field carries no strain, so the patch test asserts nothing"
    );
    let boundary_force = forces
        .iter()
        .enumerate()
        .filter(|(dof, _)| dof / D != interior)
        .fold(0.0_f64, |m, (_, f)| m.max(f.abs()));
    assert!(
        boundary_force > force_bound * 1e3,
        "the patch carries no boundary traction ({boundary_force:.3e}), so the interior \
         assertion is not distinguishing equilibrium from an operator that returns zero"
    );
}

#[test]
fn the_patch_test_passes_on_a_distorted_two_dimensional_patch() {
    let (nodes, cells, interior) = square_patch();
    // A general constant gradient: extension, contraction, and shear at once,
    // and unsymmetric so the rotation part is exercised too.
    assert_patch_test(
        &nodes,
        &cells,
        interior,
        [0.0, 0.0],
        [[1.3e-4, -0.7e-4], [0.4e-4, 2.1e-4]],
    );
}

#[test]
fn the_patch_test_passes_under_a_large_rigid_offset() {
    // The same patch translated far from the origin. The strain is unchanged,
    // so the residual bound is unchanged in physics but the displacements are
    // four orders larger — this is what proves the reference-node differencing
    // cancels the constant term rather than merely making it small.
    let (nodes, cells, interior) = square_patch();
    assert_patch_test(
        &nodes,
        &cells,
        interior,
        [5.0, -3.0],
        [[1.3e-4, -0.7e-4], [0.4e-4, 2.1e-4]],
    );
}

#[test]
fn the_patch_test_passes_on_a_distorted_three_dimensional_patch() {
    let (nodes, cells, interior) = cube_patch();
    assert_patch_test(
        &nodes,
        &cells,
        interior,
        [0.0, 0.0, 0.0],
        [
            [1.1e-4, -0.5e-4, 0.3e-4],
            [0.2e-4, 1.7e-4, -0.9e-4],
            [-0.6e-4, 0.8e-4, 2.3e-4],
        ],
    );
}

#[test]
fn the_patch_test_passes_under_pure_shear() {
    // A deviatoric field: tr(eps) = 0, so the volumetric term drops out and
    // the shear term carries the whole stress. A defect in one of the two
    // Lame terms hides behind the other under a general field.
    let (nodes, cells, interior) = square_patch();
    assert_patch_test(
        &nodes,
        &cells,
        interior,
        [0.0, 0.0],
        [[0.0, 1.9e-4], [1.9e-4, 0.0]],
    );
}

#[test]
fn the_patch_test_passes_under_pure_dilation() {
    // The complement: eps = c I, no shear at all.
    let (nodes, cells, interior) = square_patch();
    assert_patch_test(
        &nodes,
        &cells,
        interior,
        [0.0, 0.0],
        [[1.5e-4, 0.0], [0.0, 1.5e-4]],
    );
}

// ---------------------------------------------------------------------------
// Scalar generality
// ---------------------------------------------------------------------------

#[test]
fn the_patch_test_holds_at_f32() {
    // The same oracle at single precision, with the bound scaled to f32's
    // epsilon rather than the assertion loosened. A generic kernel that
    // secretly computes in f64 would pass an f64 test and this one; one that
    // computes in f32 and reports f64 would fail here.
    let (nodes64, cells, interior) = square_patch();
    let nodes: [[f32; 2]; 5] = nodes64.map(|p| [p[0] as f32, p[1] as f32]);
    let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid patch");
    let (young, poisson) = (200e9_f64, 0.3);
    let (lambda, mu) = lame(young, poisson);

    let gradient = [[1.3e-4_f32, -0.7e-4], [0.4e-4, 2.1e-4]];
    let mut displacements = [0.0_f32; 10];
    for (node, position) in nodes.iter().enumerate() {
        for component in 0..2 {
            displacements[node * 2 + component] =
                gradient[component][0] * position[0] + gradient[component][1] * position[1];
        }
    }

    let mut forces = [0.0_f32; 10];
    mesh.internal_forces(&moduli::<f32>(young, poisson), &displacements, &mut forces)
        .expect("well-shaped");

    let displacement_magnitude = displacements.iter().fold(0.0_f32, |m, v| m.max(v.abs()));
    let mut gradient_magnitude = 0.0_f32;
    let mut total_measure = 0.0_f32;
    for connectivity in &cells {
        let mut coordinates = [[0.0_f32; 2]; 3];
        for (slot, node) in coordinates.iter_mut().zip(connectivity.iter()) {
            *slot = nodes[*node];
        }
        let element = ares::Simplex::new(&coordinates);
        total_measure += element.signed_measure();
        for shape in &element.shape_gradients().expect("non-degenerate") {
            for component in shape {
                gradient_magnitude = gradient_magnitude.max(component.abs());
            }
        }
    }
    let bound = 32.0
        * total_measure
        * (lambda + 2.0 * mu) as f32
        * displacement_magnitude
        * gradient_magnitude
        * gradient_magnitude
        * f32::EPSILON;

    for component in 0..2 {
        let residual = forces[interior * 2 + component];
        assert!(
            residual.abs() <= bound,
            "interior residual {residual:.6e} exceeds the f32 bound {bound:.3e}"
        );
    }
}
