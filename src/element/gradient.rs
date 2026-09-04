use eunomia::{NumericElement, RealField};

use super::{DegenerateElement, Simplex};

impl<T: RealField, const D: usize, const N: usize> Simplex<'_, T, D, N> {
    /// The displacement gradient `grad u` this element reconstructs from its
    /// nodal displacements.
    ///
    /// Constant over the element, because the shape gradients are.
    ///
    /// # Errors
    ///
    /// Returns [`DegenerateElement`] when the element has collapsed.
    pub fn displacement_gradient(
        &self,
        displacements: &[[T; D]; N],
    ) -> Result<[[T; D]; D], DegenerateElement> {
        Ok(accumulate_gradient(&self.shape_gradients()?, displacements))
    }
}

/// `grad u = sum_a (u_a - u_0) (x) grad N_a`.
///
/// Split out so the stiffness action and the strain recovery share one
/// implementation of the differencing below — it is the load-bearing part, and
/// two copies of it would be two chances to lose it.
///
/// # Why the difference form
///
/// Mathematically identical to `sum_a u_a (x) grad N_a`, because the shape
/// gradients sum to zero. Numerically it is stronger, and that difference is
/// the whole reason for the subtraction.
///
/// The gradients cancel exactly only when summed in one order. Here they are
/// re-accumulated in another, so `sum_a grad N_a` is zero to rounding rather
/// than identically — measured at `1.4e-17` for an ordinary triangle. Under
/// the plain form a rigid translation therefore leaves a residual gradient
/// proportional to the translation, which becomes a spurious stress that grows
/// with how far the body has moved and not with the mesh, so refinement never
/// reveals it.
///
/// Differencing against node 0 removes it at the source: a uniform translation
/// makes every `u_a - u_0` exactly zero, so the gradient is exactly zero
/// whatever the shape gradients rounded to. Translation invariance becomes a
/// property of the formulation rather than of a cancellation that happens to
/// work out.
pub(crate) fn accumulate_gradient<T: RealField, const D: usize, const N: usize>(
    shape_gradients: &[[T; D]; N],
    displacements: &[[T; D]; N],
) -> [[T; D]; D] {
    let reference = displacements[0];
    let mut gradient = [[<T as NumericElement>::ZERO; D]; D];
    for (displacement, shape) in displacements.iter().zip(shape_gradients.iter()) {
        for (i, row) in gradient.iter_mut().enumerate() {
            let relative = displacement[i] - reference[i];
            for (j, entry) in row.iter_mut().enumerate() {
                *entry += relative * shape[j];
            }
        }
    }
    gradient
}
