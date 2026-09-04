//! Executable evidence for the Phase 0 constitutive coupling (atlas ADR 0057).
//!
//! Every material constant here is stated by the test, never taken from the
//! Proteus catalog. The law under test is `sigma(eps; lambda, mu)`, and
//! sourcing its inputs from a catalog would couple these assertions to
//! published values that are corrected from time to time — a catalog entry
//! moving would then fail a test of the constitutive law, which is not what
//! it measures.

#![expect(
    clippy::float_cmp,
    reason = "the closed forms below are evaluated in the same order as the implementation over exactly representable inputs, so equality is exact; the two identities with genuine rounding, the uniaxial round trip and von Mises, carry derived bounds instead."
)]

use aequitas::systems::si::quantities::{Dimensionless, Pressure};
use ares::{CauchyStress, SmallStrain, SymmetricTensor, isotropic_hooke};
use eunomia::RealField;
use proteus::IsotropicModuli;

/// Moduli from an explicit engineering pair, stated by the test.
fn moduli<T: RealField>(young: f64, poisson: f64) -> IsotropicModuli<T> {
    IsotropicModuli::from_young_poisson(
        Pressure::from_base(T::from_f64(young)),
        Dimensionless::from_base(T::from_f64(poisson)),
    )
    .expect("the stated pair is inside the positive-definite domain")
}

fn uniaxial_strain<T: RealField, const D: usize>(e: f64) -> SmallStrain<T, D> {
    let mut gradient = [[T::from_f64(0.0); D]; D];
    if let Some(row) = gradient.first_mut()
        && let Some(entry) = row.first_mut()
    {
        *entry = T::from_f64(e);
    }
    SmallStrain::from_displacement_gradient(&gradient)
}

/// Drop `//` line comments, leaving code.
///
/// Deliberately simple: the sources under test carry no string literal
/// containing `//`, so a full lexer would be machinery without a purpose here.
fn strip_comments(source: &str) -> String {
    source
        .lines()
        .map(|line| line.split_once("//").map_or(line, |(code, _)| code))
        .collect::<Vec<_>>()
        .join(
            "
",
        )
}

// ---------------------------------------------------------------------------
// The law reproduces its closed form
// ---------------------------------------------------------------------------

#[test]
fn hooke_matches_the_closed_form_for_a_uniaxial_strain() {
    // eps_xx = e, all else zero. Then tr(eps) = e and
    //   sigma_xx = lambda e + 2 mu e,  sigma_yy = sigma_zz = lambda e.
    let (young, poisson, e) = (200e9_f64, 0.3_f64, 1e-3_f64);
    let m = moduli::<f64>(young, poisson);
    let lambda = *m.lame_lambda().as_base();
    let mu = *m.shear_modulus().as_base();

    let stress = isotropic_hooke(&m, &uniaxial_strain::<f64, 3>(e));
    let c = stress.tensor().components();

    assert_eq!(c[0][0], 2.0 * mu * e + lambda * e);
    assert_eq!(c[1][1], lambda * e);
    assert_eq!(c[2][2], lambda * e);
    assert_eq!(c[0][1], 0.0);
    assert_eq!(c[1][2], 0.0);
}

#[test]
fn pure_shear_strain_gives_pure_shear_stress() {
    // eps_xy = g/2 with no diagonal: tr(eps) = 0, so sigma = 2 mu eps and the
    // hydrostatic part vanishes. This is the case a Voigt factor-of-two error
    // gets wrong.
    let m = moduli::<f64>(70e9, 0.33);
    let mu = *m.shear_modulus().as_base();
    let g = 4e-4_f64;

    // grad u with only du_x/dy = g gives eps_xy = g/2.
    let gradient = [[0.0_f64, g], [0.0, 0.0]];
    let strain = SmallStrain::<f64, 2>::from_displacement_gradient(&gradient);
    let stress = isotropic_hooke(&m, &strain);
    let c = stress.tensor().components();

    assert_eq!(c[0][0], 0.0);
    assert_eq!(c[1][1], 0.0);
    assert_eq!(c[0][1], 2.0 * mu * (g / 2.0));
    assert_eq!(c[0][1], c[1][0]);
}

#[test]
fn hydrostatic_strain_gives_hydrostatic_stress() {
    // eps = e I in 3-D: tr = 3e, so sigma = (3 lambda + 2 mu) e I = 3 K e I.
    let m = moduli::<f64>(200e9, 0.3);
    let lambda = *m.lame_lambda().as_base();
    let mu = *m.shear_modulus().as_base();
    let e = 1e-4_f64;

    let gradient = [[e, 0.0, 0.0], [0.0, e, 0.0], [0.0, 0.0, e]];
    let strain = SmallStrain::<f64, 3>::from_displacement_gradient(&gradient);
    let stress = isotropic_hooke(&m, &strain);
    let c = stress.tensor().components();

    let expected = 2.0 * mu * e + lambda * (3.0 * e);
    assert_eq!(c[0][0], expected);
    assert_eq!(c[1][1], expected);
    assert_eq!(c[2][2], expected);
    // Purely hydrostatic: no deviatoric part, so no von Mises stress.
    assert_eq!(stress.von_mises().into_base(), 0.0);
}

fn assert_rigid_motion_is_stress_free<T: RealField>() {
    // Rigid-body stress-freedom is inherited from the kinematics, not
    // re-established here: zero strain maps to zero stress by linearity.
    let w = T::from_f64(0.37);
    let zero = T::from_f64(0.0);
    let gradient = [[zero, w, -w], [-w, zero, w], [w, -w, zero]];
    let strain = SmallStrain::<T, 3>::from_displacement_gradient(&gradient);
    let stress = isotropic_hooke(&moduli::<T>(200e9, 0.3), &strain);

    for row in stress.tensor().components() {
        for entry in row {
            assert!(*entry == zero, "rigid rotation produced stress");
        }
    }
}

#[test]
fn rigid_motion_produces_exactly_zero_stress() {
    assert_rigid_motion_is_stress_free::<f32>();
    assert_rigid_motion_is_stress_free::<f64>();
}

// ---------------------------------------------------------------------------
// Independent oracle: the uniaxial round trip through Young's modulus
// ---------------------------------------------------------------------------

#[test]
fn uniaxial_stress_recovers_youngs_modulus() {
    // An independent check on the (lambda, mu) pair: under *uniaxial stress*
    // (not uniaxial strain), sigma_xx / eps_xx must equal E. Building that
    // state needs the transverse contraction eps_yy = eps_zz = -nu eps_xx,
    // which is the definition of Poisson's ratio — so this closes the loop
    // from E and nu, through Proteus's conversion, back to E.
    let (young, poisson, e) = (200e9_f64, 0.3_f64, 1e-3_f64);
    let m = moduli::<f64>(young, poisson);

    let lateral = -poisson * e;
    let gradient = [[e, 0.0, 0.0], [0.0, lateral, 0.0], [0.0, 0.0, lateral]];
    let strain = SmallStrain::<f64, 3>::from_displacement_gradient(&gradient);
    let stress = isotropic_hooke(&m, &strain);
    let c = stress.tensor().components();

    let axial = c[0][0];
    let recovered = axial / e;
    // Several multiplications and a subtraction with cancellation in the
    // transverse terms; 32 ULP at the modulus scale covers it.
    let tolerance = young * f64::EPSILON * 32.0;
    assert!(
        (recovered - young).abs() <= tolerance,
        "recovered E = {recovered}, expected {young}"
    );
    // The transverse stresses vanish, which is what makes this uniaxial
    // *stress* rather than uniaxial strain.
    assert!(c[1][1].abs() <= young * f64::EPSILON * 32.0);
    assert!(c[2][2].abs() <= young * f64::EPSILON * 32.0);
}

// ---------------------------------------------------------------------------
// Stress invariants
// ---------------------------------------------------------------------------

#[test]
fn von_mises_matches_the_uniaxial_closed_form() {
    // Under uniaxial stress the von Mises equivalent equals the axial stress.
    let (young, poisson, e) = (200e9_f64, 0.3_f64, 1e-3_f64);
    let m = moduli::<f64>(young, poisson);
    let lateral = -poisson * e;
    let gradient = [[e, 0.0, 0.0], [0.0, lateral, 0.0], [0.0, 0.0, lateral]];
    let stress = isotropic_hooke(
        &m,
        &SmallStrain::<f64, 3>::from_displacement_gradient(&gradient),
    );

    let axial = stress.tensor().components()[0][0];
    let equivalent = stress.von_mises().into_base();
    let tolerance = axial.abs() * f64::EPSILON * 64.0;
    assert!(
        (equivalent - axial).abs() <= tolerance,
        "von Mises {equivalent} differs from the axial stress {axial}"
    );
}

#[test]
fn von_mises_under_pure_shear_is_root_three_times_the_shear_stress() {
    // A second, independent witness to the shear factor. Mutation testing
    // showed only one test caught scaling shear by mu instead of 2 mu — the
    // uniaxial and hydrostatic cases carry no shear at all, exactly as the
    // module doc warns. This adds a closed form that does:
    // for pure shear tau, the von Mises equivalent is sqrt(3) tau.
    let m = moduli::<f64>(70e9, 0.33);
    let mu = *m.shear_modulus().as_base();
    let g = 4e-4_f64;

    let strain = SmallStrain::<f64, 2>::from_displacement_gradient(&[[0.0, g], [0.0, 0.0]]);
    let stress = isotropic_hooke(&m, &strain);

    let tau = stress.tensor().components()[0][1];
    assert_eq!(tau, mu * g); // 2 mu * (g/2)

    let expected = 3.0_f64.sqrt() * tau;
    let equivalent = stress.von_mises().into_base();
    let tolerance = expected.abs() * f64::EPSILON * 16.0;
    assert!(
        (equivalent - expected).abs() <= tolerance,
        "von Mises {equivalent} differs from sqrt(3) tau = {expected}"
    );
}

#[test]
fn von_mises_is_insensitive_to_hydrostatic_pressure() {
    // The defining property: adding a pressure changes no deviatoric part.
    let m = moduli::<f64>(200e9, 0.3);
    let base = isotropic_hooke(&m, &uniaxial_strain::<f64, 3>(1e-3));

    let p = 5e8_f64;
    let mut shifted = *base.tensor().components();
    for (i, row) in shifted.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            if i == j {
                *entry += p;
            }
        }
    }
    let shifted = CauchyStress::<f64, 3>::from_tensor(
        SymmetricTensor::try_from_components(shifted).expect("still symmetric"),
    );

    let before = base.von_mises().into_base();
    let after = shifted.von_mises().into_base();
    assert!(
        (after - before).abs() <= before * f64::EPSILON * 64.0,
        "hydrostatic shift changed von Mises from {before} to {after}"
    );
}

#[test]
fn the_deviatoric_stress_is_trace_free() {
    let m = moduli::<f64>(70e9, 0.33);
    let stress = isotropic_hooke(&m, &uniaxial_strain::<f64, 3>(2e-3));
    let trace = stress.deviatoric().trace();
    let scale = stress.tensor().components()[0][0].abs().max(1.0);
    assert!(trace.abs() <= scale * f64::EPSILON * 8.0);
}

#[test]
fn zero_strain_gives_zero_stress_and_zero_equivalent() {
    let m = moduli::<f64>(200e9, 0.3);
    let stress = isotropic_hooke(&m, &SmallStrain::<f64, 3>::zero());
    assert_eq!(stress.von_mises().into_base(), 0.0);
    assert_eq!(stress.mean_stress().into_base(), 0.0);
}

// ---------------------------------------------------------------------------
// Ownership: no material data lives here
// ---------------------------------------------------------------------------

#[test]
fn the_crate_names_no_material() {
    // ADR 0055 R2: a balance package stores no material constant and names no
    // alloy. Asserted against the source rather than trusted, because the
    // failure is additive — someone adds a convenience constructor for steel
    // and nothing else complains.
    //
    // Comments are stripped first. The rule governs code, not prose: the
    // module doc legitimately says that a caller wanting steel asks Proteus
    // for it, and a check that flagged its own explanation would be measuring
    // the wrong thing. The first draft of this test did exactly that.
    let code = strip_comments(include_str!("../src/constitutive/hooke.rs"));
    for name in ["steel", "aluminium", "aluminum", "titanium", "concrete"] {
        assert!(
            !code.to_lowercase().contains(name),
            "the constitutive law names {name} in code; material data belongs to Proteus"
        );
    }
}

// ---------------------------------------------------------------------------
// Properties
// ---------------------------------------------------------------------------

proptest::proptest! {
    #[test]
    fn the_law_is_linear_in_strain(
        young in 1e6_f64..5e11,
        poisson in -0.4_f64..0.45,
        e in -1e-2_f64..1e-2,
        k in -4.0_f64..4.0,
    ) {
        // sigma(k eps) = k sigma(eps): linearity is what makes the assembled
        // stiffness constant, so a nonlinearity here would silently make every
        // linear solve wrong.
        let m = moduli::<f64>(young, poisson);
        let single = isotropic_hooke(&m, &uniaxial_strain::<f64, 3>(e));
        let scaled = isotropic_hooke(&m, &uniaxial_strain::<f64, 3>(k * e));

        let expected = single.tensor().components()[0][0] * k;
        let actual = scaled.tensor().components()[0][0];
        let scale = expected.abs().max(1.0);
        proptest::prop_assert!((actual - expected).abs() <= scale * f64::EPSILON * 16.0);
    }

    #[test]
    fn stress_is_always_symmetric(
        young in 1e6_f64..5e11,
        poisson in -0.4_f64..0.45,
        a in -1e-3_f64..1e-3,
        b in -1e-3_f64..1e-3,
        d in -1e-3_f64..1e-3,
    ) {
        // Angular-momentum balance requires a symmetric Cauchy stress.
        let m = moduli::<f64>(young, poisson);
        let strain = SmallStrain::<f64, 2>::from_displacement_gradient(&[[a, b], [b, d]]);
        let stress = isotropic_hooke(&m, &strain);
        let c = stress.tensor().components();
        proptest::prop_assert_eq!(c[0][1], c[1][0]);
    }

    #[test]
    fn von_mises_is_non_negative(
        young in 1e6_f64..5e11,
        poisson in -0.4_f64..0.45,
        a in -1e-3_f64..1e-3,
        b in -1e-3_f64..1e-3,
        d in -1e-3_f64..1e-3,
    ) {
        // J2 is half a sum of squares, so its root is real and non-negative.
        // A NaN here would mean a negative J2, which is structurally
        // impossible — so this also guards the invariant's construction.
        let m = moduli::<f64>(young, poisson);
        let strain = SmallStrain::<f64, 2>::from_displacement_gradient(&[[a, b], [b, d]]);
        let equivalent = isotropic_hooke(&m, &strain).von_mises().into_base();
        proptest::prop_assert!(equivalent >= 0.0 && equivalent.is_finite());
    }
}
