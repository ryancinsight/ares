use eunomia::RealField;
use proteus::IsotropicModuli;

use super::{MisshapedField, SimplexMesh};
use crate::boundary::DirichletConditions;

impl<T: RealField, const D: usize, const N: usize> SimplexMesh<'_, T, D, N> {
    /// Apply the Dirichlet-constrained stiffness: `A = P K P + (I - P)`.
    ///
    /// `P` zeroes the constrained degrees of freedom. Writing the operator
    /// this way keeps the constrained system square over the whole field, so
    /// the solver never renumbers degrees of freedom and the solution vector
    /// is indexed the same before and after the solve.
    ///
    /// # Why this form rather than deleting rows
    ///
    /// The textbook alternative — striking out constrained rows and columns —
    /// produces a smaller system with a different numbering, so every vector
    /// crossing the boundary needs mapping in both directions. Two mappings
    /// that must agree is one more chance to disagree.
    ///
    /// This form instead leaves identity rows in place: the constrained part
    /// of the output is the constrained part of the input. `A` is then still
    /// symmetric, and still positive definite once the conditions remove the
    /// rigid-body modes, so conjugate gradients applies unchanged. The
    /// constrained block contributes the identity's eigenvalue of one, which
    /// changes the conditioning but not the solvability.
    ///
    /// Symmetry needs the projection on **both** sides. Applying it only to
    /// the output would leave the coupling column `K_cf` intact, giving a
    /// non-symmetric operator on which conjugate gradients has no convergence
    /// guarantee — and which still returns plausible numbers.
    ///
    /// # Scratch
    ///
    /// `scratch` holds the projected input, so `input` stays untouched. It is
    /// a parameter rather than an owned buffer because this crate does not
    /// allocate; a caller that applies the operator repeatedly allocates once
    /// and reuses it, which is the Krylov case.
    ///
    /// # Errors
    ///
    /// Returns [`MisshapedField`] when any buffer's length is not
    /// [`degrees_of_freedom`](SimplexMesh::degrees_of_freedom).
    pub fn constrained_action(
        &self,
        moduli: &IsotropicModuli<T>,
        conditions: &DirichletConditions<'_, T, D>,
        input: &[T],
        output: &mut [T],
        scratch: &mut [T],
    ) -> Result<(), MisshapedField> {
        let expected = self.degrees_of_freedom();
        if scratch.len() != expected {
            return Err(MisshapedField::Scratch {
                expected,
                found: scratch.len(),
            });
        }
        if input.len() != expected {
            return Err(MisshapedField::Displacement {
                expected,
                found: input.len(),
            });
        }

        scratch.copy_from_slice(input);
        conditions.project(scratch);
        self.internal_forces(moduli, scratch, output)?;
        conditions.project(output);
        conditions.carry(input, output);
        Ok(())
    }

    /// Build the right-hand side the constrained operator solves against:
    /// `b = P (f_ext - K u_g) + (I - P) g`.
    ///
    /// `external` is the applied load — traction and body force — and `g` the
    /// prescribed displacements, which the conditions carry.
    ///
    /// # Why the load carries a stiffness term
    ///
    /// A non-zero prescribed displacement does work on the free degrees of
    /// freedom through the coupling block `K_fc`. Moving it to the right-hand
    /// side is what makes the constrained system equivalent to the original
    /// one; omitting it silently solves a different problem — the one where
    /// every prescribed displacement is zero — and that problem has a
    /// perfectly convergent solution, so nothing reports the substitution.
    ///
    /// With every prescribed value zero the term vanishes and this reduces to
    /// projecting the load, which is why the omission survives so many test
    /// suites: the fixed-at-zero case is the common one.
    ///
    /// # Errors
    ///
    /// Returns [`MisshapedField`] when any buffer's length is not
    /// [`degrees_of_freedom`](SimplexMesh::degrees_of_freedom).
    pub fn constrained_load(
        &self,
        moduli: &IsotropicModuli<T>,
        conditions: &DirichletConditions<'_, T, D>,
        external: &[T],
        load: &mut [T],
        scratch: &mut [T],
    ) -> Result<(), MisshapedField> {
        let expected = self.degrees_of_freedom();
        if scratch.len() != expected {
            return Err(MisshapedField::Scratch {
                expected,
                found: scratch.len(),
            });
        }
        if external.len() != expected {
            return Err(MisshapedField::Load {
                expected,
                found: external.len(),
            });
        }

        // scratch = u_g: the prescribed values, zero elsewhere.
        for slot in scratch.iter_mut() {
            *slot = <T as eunomia::NumericElement>::ZERO;
        }
        conditions.impose(scratch);

        // load = K u_g, then b = f_ext - K u_g on the free rows.
        self.internal_forces(moduli, scratch, load)?;
        for (slot, applied) in load.iter_mut().zip(external.iter()) {
            *slot = *applied - *slot;
        }
        conditions.project(load);
        conditions.impose(load);
        Ok(())
    }
}
