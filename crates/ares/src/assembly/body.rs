use eunomia::{NumericElement, RealField};

use super::{MisshapedField, SimplexMesh};
use crate::element::Simplex;

impl<T: RealField, const D: usize, const N: usize> SimplexMesh<'_, T, D, N> {
    /// Add the consistent nodal loads for a body force sampled at the nodes.
    ///
    /// `body_force` is a force per unit volume in the flat nodal layout, and
    /// `loads` is added to rather than assigned, so a body force composes with
    /// tractions and with itself.
    ///
    /// # The consistent load
    ///
    /// The load is `f_a = integral(N_a b) dV` with `b` taken as its own linear
    /// interpolant, which for a simplex gives
    ///
    /// ```text
    /// integral(N_a N_b) dV = |Om| (1 + delta_ab) / ((D + 1)(D + 2))
    /// ```
    ///
    /// and so, summing over `b`,
    ///
    /// ```text
    /// f_a = |Om| (sum_b b_b + b_a) / ((D + 1)(D + 2))
    /// ```
    ///
    /// This is the integral the weak form actually specifies, evaluated
    /// exactly for a body force linear over the element. That is the reason
    /// for it — not an accuracy claim.
    ///
    /// # What lumping would and would not cost
    ///
    /// The lumped alternative `|Om| b_a / (D + 1)` is the vertex quadrature
    /// rule, which is **also exact for a linear integrand**, so it is second
    /// order too and does not degrade the convergence order of linear
    /// elements. An earlier version of this comment claimed it would cap
    /// refinement at `O(h)`; that was wrong, and substituting the lumped form
    /// leaves the manufactured-solution rate study at second order unchanged.
    ///
    /// Measured on that same problem at 16 divisions, the lumped load is in
    /// fact the *more* accurate of the two — 3.8e-3 relative against 1.2e-2 —
    /// because its error partly cancels the discretisation's rather than
    /// adding to it. That cancellation is a property of this problem and this
    /// mesh, not a general one, which is exactly why it is not the criterion.
    ///
    /// The consistent form is kept because it is the Galerkin statement rather
    /// than an approximation of it, and because work conservation across a
    /// coupling interface — the reason ADR 0059 will care — is a property of
    /// the consistent load and not of the lumped one.
    ///
    /// # Errors
    ///
    /// Returns [`MisshapedField`] when either buffer's length is not
    /// [`degrees_of_freedom`](SimplexMesh::degrees_of_freedom).
    pub fn add_body_force(&self, body_force: &[T], loads: &mut [T]) -> Result<(), MisshapedField> {
        let expected = self.degrees_of_freedom();
        if body_force.len() != expected {
            return Err(MisshapedField::Load {
                expected,
                found: body_force.len(),
            });
        }
        if loads.len() != expected {
            return Err(MisshapedField::Force {
                expected,
                found: loads.len(),
            });
        }

        // (D + 1)(D + 2), built by repeated addition because `T` is a real
        // field rather than something a `usize` converts into.
        let mut nodes = <T as NumericElement>::ZERO;
        for _ in 0..N {
            nodes += <T as NumericElement>::ONE;
        }
        let divisor = nodes * (nodes + <T as NumericElement>::ONE);

        let (nodal_force, _) = body_force.as_chunks::<D>();
        let (nodal_loads, _) = loads.as_chunks_mut::<D>();
        for connectivity in self.cells() {
            let coordinates = self.gather(connectivity);
            let measure = Simplex::new(&coordinates).signed_measure();

            let mut total = [<T as NumericElement>::ZERO; D];
            for node in connectivity {
                for (slot, component) in total.iter_mut().zip(nodal_force[*node].iter()) {
                    *slot += *component;
                }
            }

            for node in connectivity {
                let own = nodal_force[*node];
                for ((slot, sum), single) in nodal_loads[*node]
                    .iter_mut()
                    .zip(total.iter())
                    .zip(own.iter())
                {
                    *slot += measure * (*sum + *single) / divisor;
                }
            }
        }
        Ok(())
    }
}
