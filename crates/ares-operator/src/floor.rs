//! The rounding floor of evaluating this operator's residual.
use ares::Simplex;
use athena_leto::LetoBackendError;
use eunomia::{FloatElement, NumericElement, RealField};

use crate::ConstrainedStiffness;

impl<T: RealField, const D: usize, const N: usize> ConstrainedStiffness<'_, T, D, N> {
    /// An upper bound on the 2-norm error of evaluating `b - A u` in `T` at
    /// the displacement that solves `A u = b`.
    ///
    /// From athena#158 on, Athena's conjugate gradients accepts a solve only
    /// when the recomputed `|b - A u|_2` meets its threshold, and that
    /// evaluation errs by up to this bound at the solution, so a threshold
    /// far below it is not reliably attainable. Use it as an absolute
    /// threshold.
    ///
    /// # Derivation
    ///
    /// The action is computed per cell as `f_a = |T_e| sigma . grad N_a`,
    /// `sigma = 2 mu eps + lambda tr(eps) I`, `eps = sym(G)`,
    /// `G = sum_a (u_a - u_0) (x) grad N_a`, then summed over the cells at
    /// each node. The shape gradients, measure and moduli are fixed values
    /// computed from the mesh, so they are the operator's coefficients and
    /// carry no error into the residual. An expression built from `+`, `-`
    /// and `*` on such inputs, evaluated with unit roundoff `u`, differs from
    /// its exact value by at most `gamma_k = k u / (1 - k u)` times the same
    /// expression with every input and coefficient replaced by its absolute
    /// value and every `-` by `+`, where `k` is the largest number of
    /// operations on one path from an input to the result. That reading holds
    /// when each product has one coefficient operand; a product of two
    /// computed values adds their operation counts instead. Every product in
    /// this kernel is a computed value times a coefficient: a shape gradient
    /// (in the gradient and in the force), the exact half of the
    /// symmetrisation, `2 mu` or `lambda` (in the stress), or the measure (in
    /// the force). The source is recalled, not verified: Higham, *Accuracy
    /// and Stability of Numerical Algorithms*, 2nd ed., chapter 3.
    ///
    /// Counting the kernel's longest path per output entry:
    ///
    /// - gradient: one difference, one product, `N` sums;
    /// - symmetrisation: one sum, an exact halving;
    /// - Hooke: a trace of `D` sums, one product, one sum, then one more sum
    ///   when the stress is symmetrised;
    /// - force: one product and `D` sums for `sigma . grad N_a`, one product
    ///   by the measure;
    /// - assembly: one sum per cell at the node;
    /// - Athena's subtraction from `b`: one.
    ///
    /// That is `k = N + 2 D + 9 + c`, with `c` the most cells sharing a node.
    /// The projection and the identity rows copy values, which is exact.
    ///
    /// Write `|A||u|` for that absolute-value evaluation, so the computed
    /// residual errs by at most `gamma_k (|b| + |A||u|)` componentwise. At
    /// the exact solution `|b| = |A u| <= |A||u|` componentwise, so the floor
    /// is at most `2 gamma_k ||A||u||_2`, which is what this returns.
    ///
    /// # Limits
    ///
    /// The bound is evaluated at `displacement`, not at the solution the
    /// solve will reach. A caller that sizes it from an earlier, looser solve
    /// gets a threshold that can be too tight, which makes the solve fail to
    /// converge. The bound covers the rounding of evaluating the residual at
    /// `displacement` and no other error of the solution. The count mirrors
    /// `ares`'s `stiffness_action`, so it must change when that kernel's
    /// operation order changes.
    ///
    /// The derivation assumes no underflow or overflow and `k u < 1`. The
    /// floor is itself evaluated in `T`, so its value carries a relative
    /// rounding of order `(n + k) u` for `n` degrees of freedom, which the
    /// bound does not include.
    ///
    /// # Examples
    ///
    /// ```
    /// use aequitas::systems::si::quantities::{Dimensionless, Pressure};
    /// use ares::{DirichletConditions, PrescribedDisplacement, SimplexMesh};
    /// use ares_operator::ConstrainedStiffness;
    /// use proteus::IsotropicModuli;
    ///
    /// let nodes = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]];
    /// let cells = [[0, 1, 2], [0, 2, 3]];
    /// let held = [PrescribedDisplacement::new(0, 0, 0.0), PrescribedDisplacement::new(0, 1, 0.0)];
    /// let mesh = SimplexMesh::try_new(&nodes, &cells).expect("valid mesh");
    /// let conditions = DirichletConditions::try_new(&held, mesh.node_count()).expect("valid");
    /// let moduli = IsotropicModuli::from_young_poisson(
    ///     Pressure::from_base(200e9),
    ///     Dimensionless::from_base(0.3),
    /// )
    /// .expect("positive definite");
    /// let operator = ConstrainedStiffness::new(mesh, moduli, conditions);
    ///
    /// // A held entry is an identity row: it adds `|u_i|` to the norm that
    /// // `2 gamma_k` scales, with `k = N + 2 D + 9 + 2 = 18` as nodes 0 and 2
    /// // lie in both cells.
    /// let mut displacement = [0.0; 8];
    /// displacement[0] = -3.0;
    /// let floor = operator.residual_floor(&displacement).expect("eight degrees of freedom");
    /// let unit_roundoff = f64::EPSILON / 2.0;
    /// let gamma = 18.0 * unit_roundoff / (1.0 - 18.0 * unit_roundoff);
    /// assert_eq!(floor, 2.0 * gamma * 3.0);
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`LetoBackendError::LengthMismatch`] when `displacement` does
    /// not have one entry per degree of freedom.
    ///
    /// # Panics
    ///
    /// Does not panic. The `expect` guards the invariant
    /// [`SimplexMesh::try_new`](ares::SimplexMesh::try_new) established by
    /// computing the shape gradients of every cell.
    pub fn residual_floor(&self, displacement: &[T]) -> Result<T, LetoBackendError> {
        let mesh = self.mesh();
        let dofs = mesh.degrees_of_freedom();
        if displacement.len() != dofs {
            return Err(LetoBackendError::LengthMismatch {
                left: dofs,
                right: displacement.len(),
            });
        }
        let zero = <T as NumericElement>::ZERO;
        let two = <T as NumericElement>::ONE + <T as NumericElement>::ONE;
        // Lame's lambda is negative for a negative Poisson ratio.
        let lambda = self.moduli().lame_lambda().as_base().abs();
        let mu = self.moduli().shear_modulus().as_base().abs();

        // |P u|: the constrained entries enter the action as zero.
        let mut magnitude: Vec<T> = displacement.iter().map(|value| value.abs()).collect();
        self.conditions().project(&mut magnitude);
        let mut bound = vec![zero; dofs];
        let mut cells_at_node = vec![0_usize; mesh.node_count()];
        {
            let (nodal, _) = magnitude.as_chunks::<D>();
            let (nodal_bound, _) = bound.as_chunks_mut::<D>();
            for connectivity in mesh.cells() {
                let mut coordinates = [[zero; D]; N];
                for (slot, node) in coordinates.iter_mut().zip(connectivity) {
                    *slot = mesh.nodes()[*node];
                }
                let element = Simplex::new(&coordinates);
                let shape = element.shape_gradients().expect(
                    "invariant: try_new computed the shape gradients of every cell's coordinates",
                );
                let measure = element.signed_measure().abs();

                // |G|: node 0's difference is exactly zero, so it is skipped.
                let reference = nodal[connectivity[0]];
                let mut gradient = [[zero; D]; D];
                for (node, shape_gradient) in connectivity.iter().zip(&shape).skip(1) {
                    for (i, row) in gradient.iter_mut().enumerate() {
                        let relative = nodal[*node][i] + reference[i];
                        for (entry, component) in row.iter_mut().zip(shape_gradient) {
                            *entry += relative * component.abs();
                        }
                    }
                }

                // |sigma| = 2 mu |eps| + lambda tr|eps| I, |eps| = sym|G|.
                let mut stress = [[zero; D]; D];
                let mut trace = zero;
                for (i, row) in stress.iter_mut().enumerate() {
                    for (j, entry) in row.iter_mut().enumerate() {
                        *entry = (gradient[i][j] + gradient[j][i]) / two;
                    }
                    trace += row[i];
                }
                for (i, row) in stress.iter_mut().enumerate() {
                    for entry in row.iter_mut() {
                        *entry *= two * mu;
                    }
                    row[i] += lambda * trace;
                }

                for (node, shape_gradient) in connectivity.iter().zip(&shape) {
                    for (total, row) in nodal_bound[*node].iter_mut().zip(&stress) {
                        let mut sum = zero;
                        for (entry, component) in row.iter().zip(shape_gradient) {
                            sum += *entry * component.abs();
                        }
                        *total += measure * sum;
                    }
                    cells_at_node[*node] += 1;
                }
            }
        }
        // Constrained rows are identity rows, contributing u_i (squared below).
        self.conditions().project(&mut bound);
        self.conditions().carry(displacement, &mut bound);
        let norm = bound
            .iter()
            .fold(zero, |total, entry| total + *entry * *entry)
            .sqrt();

        let busiest = cells_at_node.into_iter().max().unwrap_or(0);
        let terms = <T as FloatElement>::from_count(N + 2 * D + 9 + busiest);
        let unit_roundoff = T::EPSILON / two;
        let gamma = terms * unit_roundoff / (<T as NumericElement>::ONE - terms * unit_roundoff);
        Ok(two * gamma * norm)
    }
}
