use eunomia::RealField;

use super::SymmetricTensor;

/// The infinitesimal (small) strain tensor.
///
/// ```text
/// eps = (grad u + grad u^T) / 2
/// ```
///
/// Strain is dimensionless — a length gradient over a length — so its
/// components carry no unit conversion. It is a distinct type from a bare
/// [`SymmetricTensor`] so a strain cannot be passed where a stress is meant;
/// the two share an algebra but not a meaning.
///
/// # Validity
///
/// The small-strain measure is the linearisation of the Green-Lagrange strain
/// about the undeformed configuration, dropping the quadratic term. It is the
/// correct measure only while displacement gradients are small, and Phase 0
/// (atlas ADR 0057) is scoped to that regime. Finite deformation is a later
/// phase with a different measure, not a tolerance change to this one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SmallStrain<T, const D: usize> {
    tensor: SymmetricTensor<T, D>,
}

impl<T: RealField, const D: usize> SmallStrain<T, D> {
    /// Build from a displacement gradient `grad u`, where
    /// `gradient[i][j] = du_i / dx_j`.
    ///
    /// # Theorem: rigid-body motion produces exactly zero strain
    ///
    /// A rigid translation has `grad u = 0`, so the result is zero trivially.
    ///
    /// An infinitesimal rigid rotation has an antisymmetric gradient,
    /// `grad u = W` with `W[j][i] = -W[i][j]`. Symmetrising gives
    /// `(W + W^T) / 2`, whose every entry is `(w + (-w)) / 2`. In IEEE-754
    /// that is exactly `+0.0` for any finite `w`, so the strain is exactly
    /// zero rather than small.
    ///
    /// This matters beyond tidiness: a strain measure that returned a tiny
    /// nonzero value under rotation would manufacture stress out of rigid
    /// motion, and the error would scale with the rotation rather than with
    /// the mesh, so no refinement study would reveal it.
    #[must_use]
    pub fn from_displacement_gradient(gradient: &[[T; D]; D]) -> Self {
        Self {
            tensor: SymmetricTensor::from_symmetrised(gradient),
        }
    }

    /// Volumetric strain `tr(eps)`, the relative volume change to first order.
    #[must_use]
    pub fn volumetric(&self) -> T {
        self.tensor.trace()
    }

    /// Deviatoric (shape-changing) strain, `eps - tr(eps)/D * I`.
    #[must_use]
    pub fn deviatoric(&self) -> SymmetricTensor<T, D> {
        self.tensor.deviator()
    }

    /// The zero strain state.
    #[must_use]
    pub fn zero() -> Self {
        Self {
            tensor: SymmetricTensor::zero(),
        }
    }
}

impl<T, const D: usize> SmallStrain<T, D> {
    /// Borrow the underlying symmetric tensor.
    #[must_use]
    pub const fn tensor(&self) -> &SymmetricTensor<T, D> {
        &self.tensor
    }
}
