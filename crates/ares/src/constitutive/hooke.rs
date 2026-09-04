use eunomia::{NumericElement, RealField};
use proteus::IsotropicModuli;

use super::CauchyStress;
use crate::kinematics::{SmallStrain, SymmetricTensor};

/// Isotropic linear-elastic stress from a small strain.
///
/// ```text
/// sigma = lambda tr(eps) I + 2 mu eps
/// ```
///
/// # Ownership
///
/// The moduli come from [`proteus::IsotropicModuli`], which owns the
/// `(E, nu) <-> (lambda, mu) <-> (c_p, c_s)` conversion contract and the named
/// material catalog. Ares stores no material constant and names no alloy
/// (atlas ADR 0055 R2): it applies a closure it is handed, and a caller that
/// wants steel asks Proteus for steel.
///
/// # Theorem
///
/// The map is linear in strain and preserves symmetry: `tr(eps)` is a scalar
/// so the first term is diagonal, and the second scales a symmetric tensor.
/// It therefore returns a valid Cauchy stress for any valid strain, with no
/// failure mode of its own — the validity that matters was established when
/// the moduli were constructed, inside the positive-definite domain.
///
/// A rigid motion has exactly zero strain
/// ([`SmallStrain::from_displacement_gradient`]), so it produces exactly zero
/// stress here: `lambda * 0 * I + 2 mu * 0`. Rigid-body stress-freedom is
/// therefore a property of the kinematics, carried through unchanged rather
/// than re-established.
#[must_use]
pub fn isotropic_hooke<T: RealField, const D: usize>(
    moduli: &IsotropicModuli<T>,
    strain: &SmallStrain<T, D>,
) -> CauchyStress<T, D> {
    let two = <T as NumericElement>::ONE + <T as NumericElement>::ONE;
    let lambda = *moduli.lame_lambda().as_base();
    let mu = *moduli.shear_modulus().as_base();
    let volumetric = strain.volumetric();

    let mut components = *strain.tensor().components();
    for (i, row) in components.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry *= two * mu;
            if i == j {
                *entry += lambda * volumetric;
            }
        }
    }

    // Scaling a symmetric tensor and adding a diagonal term preserves
    // symmetry, so this cannot fail; `from_symmetrised` is a no-op on an
    // already-symmetric array and states that invariant rather than
    // asserting it.
    CauchyStress::from_tensor(SymmetricTensor::from_symmetrised(&components))
}
