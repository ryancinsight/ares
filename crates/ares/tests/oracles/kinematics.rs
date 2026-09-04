//! Executable evidence for the Phase 0 kinematics (atlas ADR 0057).
//!
//! The rigid-body oracle is asserted **exactly**, not within a tolerance. A
//! strain measure that returned a tiny nonzero value under rigid rotation
//! would manufacture stress out of rigid motion, and the error would scale
//! with the rotation rather than the mesh — so no refinement study would
//! reveal it. Exactness is available here, so a tolerance would be strictly
//! weaker evidence.

#![expect(
    clippy::float_cmp,
    reason = "exactness is the property under test, not an accident of it: the rigid-body oracles assert that an antisymmetric gradient yields identically zero, which an epsilon comparison could not distinguish from a measure that merely happens to be small. The remaining comparisons use integers and halves that are exactly representable. The one identity with genuine rounding, the trace-free deviator, carries a derived bound instead."
)]

use ares::{SmallStrain, SymmetricTensor};
use eunomia::RealField;

// ---------------------------------------------------------------------------
// Rigid-body motion produces exactly zero strain
// ---------------------------------------------------------------------------

fn assert_translation_is_strain_free<T: RealField>() {
    // A rigid translation has a zero displacement gradient.
    let strain = SmallStrain::<T, 3>::from_displacement_gradient(&[[T::from_f64(0.0); 3]; 3]);
    for row in strain.tensor().components() {
        for entry in row {
            assert!(*entry == T::from_f64(0.0), "translation produced strain");
        }
    }
}

fn assert_infinitesimal_rotation_is_strain_free<T: RealField>(w: f64) {
    // grad u = W, antisymmetric: the infinitesimal rotation generator.
    let zero = T::from_f64(0.0);
    let omega = T::from_f64(w);
    let gradient = [
        [zero, omega, -omega],
        [-omega, zero, omega],
        [omega, -omega, zero],
    ];
    let strain = SmallStrain::<T, 3>::from_displacement_gradient(&gradient);

    for row in strain.tensor().components() {
        for entry in row {
            assert!(
                *entry == zero,
                "rotation of {w} produced nonzero strain; a rigid motion must \
                 produce exactly zero, or the solver manufactures stress from it"
            );
        }
    }
    assert!(strain.volumetric() == zero, "rotation changed volume");
}

#[test]
fn rigid_translation_produces_exactly_zero_strain() {
    assert_translation_is_strain_free::<f32>();
    assert_translation_is_strain_free::<f64>();
}

#[test]
fn infinitesimal_rotation_produces_exactly_zero_strain() {
    // Across magnitudes, including values whose halving is not exact and
    // whose square would underflow — the cancellation is exact regardless.
    for w in [1.0, 1e-8, 1e8, 0.1, 7.3, 1e-30] {
        assert_infinitesimal_rotation_is_strain_free::<f32>(w);
        assert_infinitesimal_rotation_is_strain_free::<f64>(w);
    }
}

// ---------------------------------------------------------------------------
// Strain reproduces its definition
// ---------------------------------------------------------------------------

#[test]
fn strain_is_the_symmetric_part_of_the_gradient() {
    // A general gradient: eps_ij = (g_ij + g_ji) / 2.
    let gradient = [[1.0_f64, 4.0, 6.0], [2.0, 3.0, 8.0], [0.0, 4.0, 5.0]];
    let strain = SmallStrain::<f64, 3>::from_displacement_gradient(&gradient);
    let components = strain.tensor().components();

    assert_eq!(components[0][0], 1.0);
    assert_eq!(components[1][1], 3.0);
    assert_eq!(components[2][2], 5.0);
    assert_eq!(components[0][1], 3.0); // (4 + 2) / 2
    assert_eq!(components[0][2], 3.0); // (6 + 0) / 2
    assert_eq!(components[1][2], 6.0); // (8 + 4) / 2
}

#[test]
fn strain_is_symmetric_by_construction() {
    let gradient = [[1.0_f64, 4.0, 6.0], [2.0, 3.0, 8.0], [0.0, 4.0, 5.0]];
    let strain = SmallStrain::<f64, 3>::from_displacement_gradient(&gradient);
    let c = strain.tensor().components();
    for (i, row) in c.iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            assert_eq!(*entry, c[j][i], "({i}, {j}) is not symmetric");
        }
    }
}

#[test]
fn uniaxial_extension_has_the_expected_volumetric_strain() {
    // A pure stretch along x by e: tr(eps) = e to first order.
    let e = 1e-3_f64;
    let mut gradient = [[0.0_f64; 3]; 3];
    gradient[0][0] = e;
    let strain = SmallStrain::<f64, 3>::from_displacement_gradient(&gradient);
    assert_eq!(strain.volumetric(), e);
}

// ---------------------------------------------------------------------------
// Tensor algebra
// ---------------------------------------------------------------------------

#[test]
fn the_deviator_is_trace_free() {
    // The defining property: removing the mean leaves no volumetric part.
    let gradient = [[3.0_f64, 1.0, 0.0], [1.0, 6.0, 2.0], [0.0, 2.0, 9.0]];
    let strain = SmallStrain::<f64, 3>::from_displacement_gradient(&gradient);
    let deviator = strain.deviatoric();
    // tr = 18, mean = 6, so the diagonal becomes (-3, 0, 3) and sums to zero.
    assert_eq!(deviator.trace(), 0.0);
}

#[test]
fn the_deviator_of_a_pure_dilation_vanishes() {
    // A pure volume change has no shape change.
    let e = 0.25_f64;
    let gradient = [[e, 0.0, 0.0], [0.0, e, 0.0], [0.0, 0.0, e]];
    let strain = SmallStrain::<f64, 3>::from_displacement_gradient(&gradient);
    for row in strain.deviatoric().components() {
        for entry in row {
            assert_eq!(*entry, 0.0);
        }
    }
}

#[test]
fn double_dot_pairs_every_component() {
    // A : B sums the elementwise product over both indices, off-diagonals
    // counted once per position rather than once per pair — the property a
    // Voigt representation gets wrong by a factor of two.
    let a = SmallStrain::<f64, 2>::from_displacement_gradient(&[[1.0, 2.0], [2.0, 3.0]]);
    let b = SmallStrain::<f64, 2>::from_displacement_gradient(&[[4.0, 5.0], [5.0, 6.0]]);
    // 1*4 + 2*5 + 2*5 + 3*6 = 4 + 10 + 10 + 18 = 42
    assert_eq!(a.tensor().double_dot(b.tensor()), 42.0);
}

#[test]
fn trace_sums_the_diagonal() {
    let strain = SmallStrain::<f64, 3>::from_displacement_gradient(&[
        [1.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 4.0],
    ]);
    assert_eq!(strain.tensor().trace(), 7.0);
}

// ---------------------------------------------------------------------------
// Construction boundaries
// ---------------------------------------------------------------------------

#[test]
fn an_asymmetric_array_is_rejected_rather_than_symmetrised() {
    // Silently symmetrising would let a caller believe the tensor they built
    // is the one they passed.
    let rejected = SymmetricTensor::<f64, 2>::try_from_components([[1.0, 2.0], [3.0, 4.0]])
        .expect_err("asymmetric input must be rejected");
    assert_eq!(rejected.row, 0);
    assert_eq!(rejected.column, 1);
}

#[test]
fn a_symmetric_array_is_accepted_unchanged() {
    let accepted = SymmetricTensor::<f64, 2>::try_from_components([[1.0, 2.0], [2.0, 4.0]])
        .expect("symmetric input is valid");
    assert_eq!(accepted.components()[0][1], 2.0);
    assert_eq!(accepted.trace(), 5.0);
}

#[test]
fn out_of_range_components_return_none() {
    let tensor = SymmetricTensor::<f64, 2>::zero();
    assert!(tensor.component(0, 0).is_some());
    assert!(tensor.component(2, 0).is_none());
    assert!(tensor.component(0, 2).is_none());
}

#[test]
fn the_zero_tensor_is_zero_everywhere() {
    let tensor = SymmetricTensor::<f64, 3>::zero();
    assert_eq!(tensor.trace(), 0.0);
    assert_eq!(tensor.double_dot(&tensor), 0.0);
}

// ---------------------------------------------------------------------------
// Properties
// ---------------------------------------------------------------------------

proptest::proptest! {
    #[test]
    fn any_antisymmetric_gradient_gives_exactly_zero_strain(w in -1e6_f64..1e6) {
        let gradient = [[0.0, w], [-w, 0.0]];
        let strain = SmallStrain::<f64, 2>::from_displacement_gradient(&gradient);
        for row in strain.tensor().components() {
            for entry in row {
                proptest::prop_assert_eq!(*entry, 0.0);
            }
        }
    }

    #[test]
    fn strain_is_always_symmetric(
        a in -1e3_f64..1e3, b in -1e3_f64..1e3,
        c in -1e3_f64..1e3, d in -1e3_f64..1e3,
    ) {
        let strain = SmallStrain::<f64, 2>::from_displacement_gradient(&[[a, b], [c, d]]);
        let components = strain.tensor().components();
        proptest::prop_assert_eq!(components[0][1], components[1][0]);
    }

    #[test]
    fn the_deviator_is_always_trace_free(
        a in -1e3_f64..1e3, b in -1e3_f64..1e3, d in -1e3_f64..1e3,
    ) {
        let strain = SmallStrain::<f64, 2>::from_displacement_gradient(&[[a, b], [b, d]]);
        let trace = strain.deviatoric().trace();
        // Removing the mean is two subtractions, so the residual is bounded by
        // the rounding of the trace itself rather than being exactly zero.
        let scale = a.abs().max(d.abs()).max(1.0);
        proptest::prop_assert!(trace.abs() <= f64::EPSILON * 8.0 * scale);
    }

    #[test]
    fn double_dot_is_commutative(
        a in -1e3_f64..1e3, b in -1e3_f64..1e3, d in -1e3_f64..1e3,
        e in -1e3_f64..1e3, f in -1e3_f64..1e3, g in -1e3_f64..1e3,
    ) {
        let x = SmallStrain::<f64, 2>::from_displacement_gradient(&[[a, b], [b, d]]);
        let y = SmallStrain::<f64, 2>::from_displacement_gradient(&[[e, f], [f, g]]);
        proptest::prop_assert_eq!(
            x.tensor().double_dot(y.tensor()),
            y.tensor().double_dot(x.tensor())
        );
    }
}
