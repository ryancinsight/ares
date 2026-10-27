use aequitas::Quantity;
use aequitas::dimension::MultiplyDimension;
use aequitas::systems::si::dimensions::Pressure as PressureDimension;
use aequitas::systems::si::quantities::Pressure;
use eunomia::{NumericElement, RealField};

use crate::kinematics::SymmetricTensor;

/// Second-deviatoric-invariant quantity, in pressure squared.
///
/// Composed through [`MultiplyDimension`] rather than declared by hand, so the
/// dimension algebra stays the single source: `PressureSquared` is
/// `Pressure × Pressure`, and `PressureSquared::sqrt()` is `Pressure` by the
/// type-level root, which [`CauchyStress::von_mises`] relies on.
pub type PressureSquared<T = f64> =
    Quantity<T, <PressureDimension as MultiplyDimension<PressureDimension>>::Output>;

/// The Cauchy stress tensor, in pascals.
///
/// A distinct type from [`SmallStrain`](crate::SmallStrain) even though both
/// wrap a symmetric tensor: they share an algebra and not a meaning, and the
/// constitutive law maps one to the other. Passing a strain where a stress is
/// meant is a type error rather than a silent factor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CauchyStress<T, const D: usize> {
    tensor: SymmetricTensor<T, D>,
}

impl<T: RealField, const D: usize> CauchyStress<T, D> {
    /// Wrap a symmetric tensor whose components are stresses in pascals.
    #[must_use]
    pub const fn from_tensor(tensor: SymmetricTensor<T, D>) -> Self {
        Self { tensor }
    }

    /// Hydrostatic (mean) stress `tr(sigma) / D`.
    ///
    /// Positive in tension, following the sign convention of the stress
    /// tensor itself. A pressure is the negative of this.
    #[must_use]
    pub fn mean_stress(&self) -> Pressure<T> {
        Pressure::from_base(self.tensor.trace() / Self::dimension())
    }

    /// Deviatoric stress `s = sigma - tr(sigma)/D * I`.
    ///
    /// The shape-changing part, which yielding depends on and hydrostatic
    /// loading does not affect.
    #[must_use]
    pub fn deviatoric(&self) -> SymmetricTensor<T, D> {
        self.tensor.deviator()
    }

    /// The second deviatoric invariant `J2 = (s : s) / 2`, in pressure squared.
    ///
    /// Closed forms: `sigma^2 / 3` under uniaxial stress `sigma`, and `tau^2`
    /// under pure shear `tau`. Hydrostatic loading contributes nothing, since
    /// the deviator of a spherical tensor vanishes.
    #[must_use]
    pub fn second_deviatoric_invariant(&self) -> PressureSquared<T> {
        let two = <T as NumericElement>::ONE + <T as NumericElement>::ONE;
        let deviator = self.deviatoric();
        PressureSquared::from_base(deviator.double_dot(&deviator) / two)
    }

    /// Von Mises equivalent stress `sqrt(3 J2)`.
    ///
    /// The product `3 * J2` stays a pressure-squared quantity and the root
    /// halves the dimension at the type level, so a change to the invariant's
    /// dimension that broke this relation would fail to compile.
    ///
    /// # Why via `J2` rather than principal stresses
    ///
    /// The equivalent stress is often written from principal stresses, which
    /// requires an eigenvalue decomposition. `sqrt(3 J2)` is the same quantity
    /// computed from invariants alone: no eigen-solve, no ordering convention,
    /// and no iteration to converge. `J2` is non-negative by construction — it
    /// is half a sum of squares — so the root is always real.
    #[must_use]
    pub fn von_mises(&self) -> Pressure<T> {
        let three =
            <T as NumericElement>::ONE + <T as NumericElement>::ONE + <T as NumericElement>::ONE;
        (self.second_deviatoric_invariant() * three).sqrt()
    }

    #[inline]
    fn dimension() -> T {
        let mut count = <T as NumericElement>::ZERO;
        for _ in 0..D {
            count += <T as NumericElement>::ONE;
        }
        count
    }
}

impl<T, const D: usize> CauchyStress<T, D> {
    /// Borrow the underlying symmetric tensor.
    #[must_use]
    pub const fn tensor(&self) -> &SymmetricTensor<T, D> {
        &self.tensor
    }
}
