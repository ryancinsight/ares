use aequitas::systems::si::quantities::Pressure;
use eunomia::{NumericElement, RealField};

use crate::kinematics::SymmetricTensor;

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

    /// The second deviatoric invariant `J2 = (s : s) / 2`.
    #[must_use]
    pub fn second_deviatoric_invariant(&self) -> T {
        let two = <T as NumericElement>::ONE + <T as NumericElement>::ONE;
        let deviator = self.deviatoric();
        deviator.double_dot(&deviator) / two
    }

    /// Von Mises equivalent stress `sqrt(3 J2)`.
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
        Pressure::from_base((three * self.second_deviatoric_invariant()).sqrt())
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
