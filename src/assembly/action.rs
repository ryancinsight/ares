use eunomia::{NumericElement, RealField};
use proteus::IsotropicModuli;

use super::SimplexMesh;
use crate::element::{Simplex, stiffness_action};
use crate::kinematics::SmallStrain;

impl<T: RealField, const D: usize, const N: usize> SimplexMesh<'_, T, D, N> {
    /// Assemble the internal force vector `f = K u` over the whole mesh.
    ///
    /// `displacements` and `forces` are flat, `D` components per node, in node
    /// order: degree of freedom `node * D + component`. `forces` is written,
    /// not accumulated into — this is an assignment, so any prior contents are
    /// overwritten.
    ///
    /// # Why a force vector rather than a stiffness matrix
    ///
    /// Athena's solvers take a `LinearOperator` and ask it for `K u`, never for
    /// `K`. Assembling a global sparse matrix would mean building a structure
    /// nothing reads: the sparsity pattern, the coordinate-to-compressed
    /// conversion, and the storage all exist only to be multiplied by a vector
    /// once per iteration, which this does directly from the mesh.
    ///
    /// It is also the residual's internal half. Static equilibrium is
    /// `f_int(u) = f_ext`, so `A6`'s solve and the balance statement in ADR
    /// 0057 are the same quantity read two ways.
    ///
    /// # Allocation
    ///
    /// None. Each cell gathers into stack buffers sized by `N`, so the loop
    /// touches the heap zero times regardless of mesh size.
    ///
    /// # Errors
    ///
    /// Returns [`MisshapedField`] when either buffer's length is not
    /// [`degrees_of_freedom`](SimplexMesh::degrees_of_freedom). No other
    /// failure is reachable: [`SimplexMesh::try_new`] already established that
    /// every cell integrates.
    ///
    /// # Panics
    ///
    /// Does not panic. The `expect` guards an invariant `try_new` established
    /// by running the identical call on the identical coordinates.
    pub fn internal_forces(
        &self,
        moduli: &IsotropicModuli<T>,
        displacements: &[T],
        forces: &mut [T],
    ) -> Result<(), MisshapedField> {
        let expected = self.degrees_of_freedom();
        if displacements.len() != expected {
            return Err(MisshapedField::Displacement {
                expected,
                found: displacements.len(),
            });
        }
        if forces.len() != expected {
            return Err(MisshapedField::Force {
                expected,
                found: forces.len(),
            });
        }

        // Both lengths are `node_count * D`, so each split leaves no remainder
        // and the nodal view costs nothing.
        let (nodal_displacements, _) = displacements.as_chunks::<D>();
        let (nodal_forces, _) = forces.as_chunks_mut::<D>();
        for slot in nodal_forces.iter_mut().flat_map(|f| f.iter_mut()) {
            *slot = <T as NumericElement>::ZERO;
        }

        let mut cell_displacements = [[<T as NumericElement>::ZERO; D]; N];
        let mut cell_forces = [[<T as NumericElement>::ZERO; D]; N];
        for connectivity in self.cells() {
            let coordinates = self.gather(connectivity);
            for (slot, node) in cell_displacements.iter_mut().zip(connectivity.iter()) {
                *slot = nodal_displacements[*node];
            }

            stiffness_action(
                &Simplex::new(&coordinates),
                moduli,
                &cell_displacements,
                &mut cell_forces,
            )
            .expect(
                "invariant: try_new ran shape_gradients on every cell's coordinates and kept only \
                 the meshes where all of them succeeded",
            );

            for (contribution, node) in cell_forces.iter().zip(connectivity.iter()) {
                for (total, part) in nodal_forces[*node].iter_mut().zip(contribution.iter()) {
                    *total += *part;
                }
            }
        }
        Ok(())
    }

    /// The constant strain each cell reconstructs from a displacement field,
    /// in cell order.
    ///
    /// Linear simplices carry constant strain, so one tensor per cell is the
    /// complete answer rather than a sampled approximation of one. Stress
    /// follows by passing each strain through
    /// [`isotropic_hooke`](crate::isotropic_hooke) — this deliberately stops
    /// at the kinematic quantity, because the constitutive step is Proteus's
    /// to close (atlas ADR 0055) and duplicating it here would give the crate
    /// two closures to keep agreeing.
    ///
    /// # Errors
    ///
    /// Returns [`MisshapedField`] when the field length is not
    /// [`degrees_of_freedom`](SimplexMesh::degrees_of_freedom).
    ///
    /// # Panics
    ///
    /// Does not panic, for the reason given on
    /// [`internal_forces`](SimplexMesh::internal_forces).
    pub fn cell_strains<'field>(
        &'field self,
        displacements: &'field [T],
    ) -> Result<impl Iterator<Item = SmallStrain<T, D>> + 'field, MisshapedField> {
        let expected = self.degrees_of_freedom();
        if displacements.len() != expected {
            return Err(MisshapedField::Displacement {
                expected,
                found: displacements.len(),
            });
        }
        let (nodal, _) = displacements.as_chunks::<D>();
        Ok(self.cells().iter().map(move |connectivity| {
            let coordinates = self.gather(connectivity);
            let mut cell_displacements = [[<T as NumericElement>::ZERO; D]; N];
            for (slot, node) in cell_displacements.iter_mut().zip(connectivity.iter()) {
                *slot = nodal[*node];
            }
            let gradient = Simplex::new(&coordinates)
                .displacement_gradient(&cell_displacements)
                .expect(
                    "invariant: try_new ran shape_gradients on every cell's coordinates and kept \
                     only the meshes where all of them succeeded",
                );
            SmallStrain::from_displacement_gradient(&gradient)
        }))
    }
}

/// A field whose length does not match the mesh.
///
/// The only failure assembly has left, and it is a caller error rather than a
/// property of the mesh.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum MisshapedField {
    /// The displacement field is the wrong length.
    Displacement {
        /// Degrees of freedom the mesh has.
        expected: usize,
        /// Length supplied.
        found: usize,
    },
    /// The force field is the wrong length.
    Force {
        /// Degrees of freedom the mesh has.
        expected: usize,
        /// Length supplied.
        found: usize,
    },
    /// The applied-load field is the wrong length.
    Load {
        /// Degrees of freedom the mesh has.
        expected: usize,
        /// Length supplied.
        found: usize,
    },
    /// The scratch buffer is the wrong length.
    Scratch {
        /// Degrees of freedom the mesh has.
        expected: usize,
        /// Length supplied.
        found: usize,
    },
}

impl core::fmt::Display for MisshapedField {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let (role, expected, found) = match self {
            Self::Displacement { expected, found } => ("displacement", expected, found),
            Self::Force { expected, found } => ("force", expected, found),
            Self::Load { expected, found } => ("applied load", expected, found),
            Self::Scratch { expected, found } => ("scratch", expected, found),
        };
        write!(
            formatter,
            "the {role} field has {found} entries, but the mesh has {expected} degrees of freedom"
        )
    }
}

impl core::error::Error for MisshapedField {}
